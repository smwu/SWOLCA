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

% sim_confounder takes in sim_n as a command line array index 
% argument, and simulates a population dataset for scenario 2.
% Inputs:
%   sim_n: Command line argument indicating array index
% Outputs: saves simulated dataset for the given scenario and index.
function sim_confounder(sim_n)
    rng(sim_n, 'twister');  % Set seed
    out_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Toy_Example/"; % Output directory
    
    %% Set scenario specifications  
    p = 4;                                      % Number of food items
    d = 2;                                      % Number of response levels (assumed constant across items)
    S = 2;                                      % Number of subpops
    K = 2;                                      % Number of global classes
    N_s = [1000, 7000];                         % Subpopulation sizes
    N = sum(N_s);                               % Population size    
    clust_mode = 0.9;                          % Probability of true consumption level occurring
    non_mode = 0.1;                            % Probability of other consumption levels occurring
       
    %% Set parameter specifications
    sim_data = struct;  % Initialize structural array
    sim_data.true_pi = [0.8, 0.2;
                        0.2, 0.8];                % Global class membership proportions
    sim_data.true_xi = [1.2, -0.33, 0.08, -0.68]; % True probit model coefficients; first dem, then class

    %% Set true comsumption patterns for each diet profile class      
    global1 = ones(p, 1) * 2;  % Global profile patterns
    global2 = ones(p, 1);
    % p x K matrix of true consumption levels for each global class
    sim_data.true_global_patterns = [global1 global2];  
    
    %% Generate true global class assignments and subpopulation assignments for all individuals
    Ci_pop = cell(S, 1);                 % Initialize global class assignments for all indivs
    Si_pop = cell(S, 1);                 % Initialize subpop assignments for all indivs
    for s = 1:S                          % For each subpopulation
        % Randomly generate global class assigns for all indivs in subpop s
        Ci_pop{s} = randsample(1:K, N_s(s), true, sim_data.true_pi(s, :))';
        Si_pop{s} = ones(N_s(s), 1) * s; % Assign subpop for all indivs in the subpop
    end    
    Si_pop = vertcat(Si_pop{:});         % Convert cell array to vector through vertical concat
    Ci_pop = vertcat(Ci_pop{:});
        
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
    q_dem = zeros(N, S);             % Initialize design matrix of dem covs (subpop) in cell-means format
    for s = 1:S                      % For each subpop
        q_dem(Si_pop == s, s) = 1;   % For all indivs in subpop s, set s-th column to 1
    end
    q_class = zeros(N, K);           % Initialize design matrix for global class membership in cell-means format
    for k = 1:K                      % For each global class
        q_class(Ci_pop == k, k) = 1; % For all indivs in global class k, set k-th column to 1
    end
    Q = [q_dem q_class];                            % Design matrix with dem covs and global classes
    lin_pred_pop = Q * transpose(sim_data.true_xi); % True linear predictor for all indivs. Mean of truncated normal dist
    Phi_pop = normcdf(lin_pred_pop);                % True probit mean, P(Y_i=1|Q, C), for all indivs
    Y_pop = binornd(1, Phi_pop);                    % True outcome for all indivs
    
    %% Format and save data
    scen = 2;  % Simulation scenario  
    n_s = N;   % Sample size is full population
    % Obtain indices, weights, and data for sampled individuals
    Li_pop = zeros(N, 1);
    sim_data = sample_indivs(N, n_s, S, false, Si_pop, Ci_pop, Li_pop, X_pop, Y_pop, Phi_pop, K, sim_data);
    % Save simulated data
    save(strcat(out_dir, 'simdata_scen', num2str(scen),'_iter', num2str(sim_n)), 'sim_data');
    
end




%% Testing code
% out_dir = "C:/Users/Lang/Documents/Harvard/Research/Briana/supRPC/wsOFMM/Toy_Example/";  % Output directory
% N_s = [200, 200, 200, 200]; 
% n_s = 40;
% n_s = [10, 10, 10, 10];
