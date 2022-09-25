%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulated data with only global latent classes and class membership 
% probabilities dependent on subpopulation 
% Programmer: SW             
% 
% Scenario 2:
% We assume individuals come from 2 subpopulations of unequal sizes. 
% There are 2 global latent classes, with class membership probabilities 
% dependent on subpopulation membership. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%54%%%%%%%%%%%%%%

% sim_confounder takes in sim_n as a command line array index argument, and 
% simulates a population dataset for that iteration. Generally, only 1 
% iteration is necessary, and samples can be generated from that dataset.
% Inputs:
%   sim_n: Command line argument indicating array index
%   scen: Command line integer indicating name for scenario number
% Outputs: saves simulated dataset for the given scenario and index.
function sim_confounder(sim_n, scen)
    rng(sim_n, 'twister');  % Set seed
    out_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Data/"; % Output directory
%             out_dir = strcat(pwd, "/");
    
    %% Set scenario specifications  
    p = 30;                                      % Number of food items
    d = 4;                                      % Number of response levels (assumed constant across items)
    S = 2;                                      % Number of subpops
    K = 3;                                      % Number of global classes
    N_s = [10000, 70000];                         % Subpopulation sizes
    N = sum(N_s);                               % Population size    
    clust_mode = 0.85;                          % Probability of true consumption level occurring
    non_mode = 0.05;                            % Probability of other consumption levels occurring
       
    %% Set parameter specifications
    sim_data = struct;  % Initialize structural array
    sim_data.true_pi = [0.2, 0.4, 0.4;
                        0.7, 0.2, 0.1];                % Global class membership proportions
    sim_data.true_xi = [1, 0, -1, 0.5, -0.5, -1.5];   % Coefficients for factor variable coding

    %% Set true comsumption patterns for each diet profile class      
    global1 = [ones(0.5 * p, 1) * 1;  % Global profile patterns
               ones(0.5 * p, 1) * 3];
    global2 = [ones(0.2 * p, 1) * 4; 
               ones(0.8 * p, 1) * 2];
    if scen == 1   % Baseline
        global3 = [ones(0.3 * p, 1) * 3; 
                   ones(0.4 * p, 1) * 2; 
                   ones(0.3 * p, 1) * 1];
    elseif scen ==2   % Supervised: global3 similar to global1
        global3 = [ones(0.1 * p, 1) * 4; 
                   ones(0.4 * p, 1) * 1; 
                   ones(0.5 * p, 1) * 3];
    end           
    
    % p x K matrix of true consumption levels for each global class
    sim_data.true_global_patterns = [global1 global2 global3];  
    
    %% Generate true global class assignments and subpopulation assignments for all individuals
    Ci_pop = cell(S, 1);                 % Initialize global class assignments for all indivs
    Si_pop = cell(S, 1);                 % Initialize subpop assignments for all indivs
    for s = 1:S                          % For each subpopulation
        pop = 1:N_s(s);                  % Initialize subpop indices
        Si_pop{s} = ones(N_s(s), 1) * s; % Assign subpop for all indivs in the subpop
        Ci_pop{s} = ones(N_s(s), 1);     % Assign base class assignment
        for c = 2:K                      % For all remaining classes
            selected = randsample(pop, N_s(s) * sim_data.true_pi(s, c)); 
            Ci_pop{s}(selected) = c;     % Assign class to sampled indivs
            pop = setdiff(pop, selected); % Redefine remaining pop to be sampled
        end
    end
    Si_pop = vertcat(Si_pop{:});         % Convert cell array to vector through vertical concat
    Ci_pop = vertcat(Ci_pop{:});
    
%             for s = 1:S                          % For each subpopulation
%                 % Randomly generate global class assigns for all indivs in subpop s
%                 Ci_pop{s} = randsample(1:K, N_s(s), true, sim_data.true_pi(s, :))';
%                 Si_pop{s} = ones(N_s(s), 1) * s; % Assign subpop for all indivs in the subpop
%             end    
%             Si_pop = vertcat(Si_pop{:});         % Convert cell array to vector through vertical concat
%             Ci_pop = vertcat(Ci_pop{:});

    %% Create true item response probabilities
    sim_data = create_item_response_probs_no_local(p, d, clust_mode, non_mode, sim_data, K);
    
    %% Create population consumption data
    X_pop = zeros(N, p);  % Initialize population consumption data for each individual and item
    for i = 1:N
        for j = 1:p
            % Create population consumption data for each indiv and item
            X_pop = create_consump_data_no_local(i, j, Ci_pop(i), sim_data, X_pop);
        end
    end  
    
    %% Create true probit model and population outcome data 
    % Factor variable coding
    q_dem = dummyvar(Si_pop);
    x_ci = dummyvar(Ci_pop);
    Q = zeros(N, S * K);          % Initialize Q design matrix
    Y_pop = zeros(N, 1);          % Initialize outcome vector
    Phi = normcdf(sim_data.true_xi);
    for s = 1:S
        for c = 1:K
            % For each (S,C) subset, use P(Y=1|S,C) to sample indivs w/ event
            subset = find((Si_pop == s) & (Ci_pop == c)); 
            with_event = randsample(subset, round(length(subset) * Phi((s-1)*K + c)));
            Y_pop(with_event) = 1;
            % Populate design matrix
            Q(:, ((s-1) * K) + c) = q_dem(:, s) .* x_ci(:, c);
        end
    end
    lin_pred_pop = Q * transpose(sim_data.true_xi); % True linear predictor for all indivs. Mean of truncated normal dist
    Phi_pop = normcdf(lin_pred_pop);                % True probit mean, P(Y_i=1|Q, C), for all indivs
%             Y_pop = binornd(1, Phi_pop);                    % True outcome for all indivs
    
    %% Format and save data
    n_s = N;   % Sample size is full population
    % Obtain indices, weights, and data for sampled individuals
    Li_pop = zeros(N, 1);
    sim_data = sample_indivs(N, n_s, S, false, Si_pop, Ci_pop, Li_pop, X_pop, Y_pop, Phi_pop, K, sim_data);
    
    disp('P(C=k)');
    tabulate(Ci_pop);
    disp('P(C=k|S)');
    tabulate(Ci_pop(1:N_s(1)));
    tabulate(Ci_pop((N_s(1)+1) : N));
    disp('P(Y=1|S,C)');
    tabulate(Phi_pop);
    disp('P(Y=1)');
    tabulate(Y_pop);

    % Save simulated data
    save(strcat(out_dir, 'simdata_scen', num2str(scen),'_iter', num2str(sim_n)), 'sim_data', 'X_pop', 'Y_pop');
    
end




%% Testing code
% % Run in command window
% sim_confounder(1, 1);  
% load(strcat(pwd, '/simdata_scen', num2str(1), '_iter', num2str(1), '.mat'), 'sim_data')

% out_dir = "C:/Users/Lang/Documents/Harvard/Research/Briana/supRPC/wsOFMM/Toy_Example/";  % Output directory
% N_s = [200, 200, 200, 200]; 
% n_s = 40;
% n_s = [10, 10, 10, 10];

% true_xi = [1.2, -0.33, 0.08, -0.68];
% Q_test = [1 0 1 0;   % S=1,C=1
%           0 1 1 0;   % S=2,C=1
%           1 0 0 1;   % S=1,C=2
%           0 1 0 1];  % S=2,C=2
% disp(Q_test * transpose(true_xi));
% disp(normcdf(Q_test * transpose(true_xi)));  
% norminv([0.8997 0.4013 0.6985 0.1562])
% 
% true_xi_ref = [1.4, -1.4, -0.5, 0.238];
% Q_test_ref = [1 0 0 0;   % S=1,C=1
%               1 1 0 0;   % S=2,C=1
%               1 0 1 0;   % S=1,C=2
%               1 1 1 1];  % S=2,C=2
% disp(Q_test_ref * transpose(true_xi_ref));
% disp(normcdf(Q_test_ref * transpose(true_xi_ref)));  
% inv(Q_test_ref) * transpose(norminv([0.9192, 0.5, 0.8159, 0.3967]))
% % Get coefs given baseline S=2,C=2:
% inv(Q_test_ref) * transpose([-0.262, 0.9, 0, 1.4])
% % Equivalently, 
% inv(Q_test_ref) * transpose(norminv([0.3967, 0.8159, 0.5, 0.9192]))

% % Get coefs given no interaction and baseline S=2,C=2
% Q_test_ref = [1 0 0;   % S=1,C=1
%               1 1 0;   % S=2,C=1
%               1 0 1;   % S=1,C=2
%               1 1 1];  % S=2,C=2
% inv(Q_test_ref) * transpose(norminv([0.3967, 0.8159, 0.5, 0.9192]))
