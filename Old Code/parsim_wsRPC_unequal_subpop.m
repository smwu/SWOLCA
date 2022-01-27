%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parallelized simulated data for weighted supervised OFMM and RPC        
% Programmer: SW             
% 
% Assume unequal subpopulation sizes
% Scenario 5: Full population. All weights equal to 1  
% Scenario 6: Sample 5% of total population (SRS). All weights equal
% Scenario 7: Sample 5% from each subpop (proportional allocation). 
%             All weights equal up to rounding      
% Scenario 8: Sample 1000 from each subpop (equal allocation). 
%             Diff weights per subpop 
%
% Data description:
% We assume individuals come from 4 subpopulations. Sampling of 4000 
% subjects is stratified by subpop. There are 3 global dietary patterns, 
% and for each subpopulation, there are 2 local dietary patterns. 
% Each food item for an individual belongs to a global pattern or to a 
% subpopulation-specific local pattern.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%54%%%%%%%%%%%%%%

clear; clc;
%% Set scenario specifications
num_sims = 100;                             % Number of simulation sets      
p = 50;                                     % Number of food items
d = 4;                                      % Number of response levels (assumed constant across items)
S = 4;                                      % Number of subpops
K = 3;                                      % Number of global classes
L = 2;                                      % Number of local classes per subpop
N_s = [10000, 25000, 27500, 17500];         % Subpopulation sizes           
N = sum(N_s);                               % Population size
clust_mode = 0.85;                          % Probability of true consumption level occurring
non_mode = 0.05;                            % Probability of other consumption levels occurring

%% Create and save simulations
parfor sim_n = 1:num_sims  % For each simulation iteration    
    rng(sim_n)     % Set seed
       
    %% Set parameter specifications
    sim_data = struct;  % Initialize structural array
    sim_data.true_pi = [0.33 0.33 0.34];        % Global class membership proportions
    sim_data.true_lambda = [0.5 0.5];           % Local class membership proportions
    sim_data.true_xi = [1 1.1 1.2 0.9 1 -2 -1]; % True probit model coefficients
    sim_data.nu = [0.5 0.4 0.7 0.8];            % Subpop-specific prob of global assigment for all items

    %% Set true comsumption patterns for each diet profile class      
    global1 = [ones(0.5 * p, 1) * 3;  % Global profile patterns
               ones(0.5 * p, 1) * 1];
    global2 = [ones(0.2 * p, 1) * 2; 
               ones(0.8 * p, 1) * 4];
    global3 = [ones(0.2 * p, 1) * 1; 
               ones(0.4 * p, 1) * 2; 
               ones(0.4 * p, 1) * 3];
    % p x K matrix of true consumption levels for each global class
    sim_data.true_global_patterns = [global1 global2 global3];  

    local11 = ones(p, 1) * 1;  % Local profile patterns
    local12 = [repmat([1; 2; 3; 4], [12, 1]); 1; 1];
    local21 = ones(p, 1) * 2;
    local22 = [repmat([2; 4; 1; 3], [12, 1]); 2; 2];
    local31 = ones(p, 1) * 3;
    local32 = [repmat([3; 1; 4; 2], [12, 1]); 3; 3];
    local41 = ones(p, 1) * 4;
    local42 = [repmat([4; 3; 2; 1], [12, 1]); 4; 4];
    % p x S x L array of true consumption levels for each local class
    sim_data.true_local_patterns = [local11 local21 local31 local41];
    sim_data.true_local_patterns(:,:,2) = [local12 local22 local32 local42];
    
    %% Generate true global and local class assignments for all individuals
    Ci_pop = randsample(1:K, N, true, sim_data.true_pi)';     % Randomly generate global class assigns for all indivs
    Li_pop = randsample(1:L, N, true, sim_data.true_lambda)'; % Randomly generate local class assigns for all indivs
    
    %% Generate subpopulation assignments for all individuals
    Si_pop = cell(S, 1);                 % Initialize subpop assignments for all indivs
    for s = 1:S                          % For each subpopulation
        Si_pop{s} = ones(N_s(s), 1) * s; % Assign subpop for all indivs in the subpop
    end    
    Si_pop = vertcat(Si_pop{:});         % Convert cell array to vector through vertical concat
        
    %% Create true item response probabilities
    sim_data = create_item_response_probs(p, d, clust_mode, non_mode, sim_data, S, K, L);
    
    %% Create deviation indicator for global or local assignment
    % Assume all indivs in the same subpop have the same items that are global/local, i.e., G_ij = G_sj
    % Assume all items for those in the same subpop have the same deviation prob, i.e., nu_sj = nu_s
    sim_data.Gsj = zeros(S, p);    % Initialize matrix of true global/local assign for each subpop and item
    for s = 1:S                    % For each subpop
        % Randomly sample vector of global/local assigns from Ber(nu(s))
        sim_data.Gsj(s, :) = binornd(1, sim_data.nu(s), [1, p]);  
    end
    nu_sum = sum(sim_data.Gsj, 2); % Number of global items for each subpop
    
    %% Create population consumption data
    X_pop = zeros(N, p);  % Initialize population consumption data for each individual and item
    for i = 1:N
        for j = 1:p
            % Create population consumption data for each indiv and item
            X_pop = create_consump_data(i, j, Si_pop(i), Ci_pop(i), Li_pop(i), sim_data, X_pop);
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
    Q = [q_dem q_class];             % Design matrix with dem covs and global classes
    lin_pred_pop = Q * transpose(sim_data.true_xi); % True linear predictor for all indivs. Mean of truncated normal dist
    sim_data.true_Phi = normcdf(lin_pred_pop);      % True probit mean, P(Y_i=1|Q, C)
    Y_pop = binornd(1, sim_data.true_Phi);          % True outcome for all indivs
    
    %% Sampling scenario 5: full population
    scen = 5;  
    n_s = N;  % Sample size
    % Obtain indices, weights, and data for sampled individuals
    sim_data = sample_indivs(N, n_s, S, false, Si_pop, Ci_pop, Li_pop, X_pop, Y_pop, K, sim_data);
    % Save simulated data
    par_save(strcat('parsimdata_wsRPC_scen', num2str(scen),'_iter', num2str(sim_n)), sim_data);
    
    %% Sampling scenario 6: simple random sample
    scen = 6;                    
    n_s = 4000;  % Sample size 
    % Obtain indices, weights, and data for sampled individuals
    sim_data = sample_indivs(N, n_s, S, false, Si_pop, Ci_pop, Li_pop, X_pop, Y_pop, K, sim_data);
    % Save simulated data
    par_save(strcat('parsimdata_wsRPC_scen', num2str(scen),'_iter', num2str(sim_n)), sim_data);
    
    %% Sampling scenario 7: stratified random sample by subpop with proportional allocation
    scen = 7;                     
    n_s = round(0.05 .* N_s);  % Sample sizes for each subpop
    % Obtain indices, weights, and data for sampled individuals
    sim_data = sample_indivs(N_s, n_s, S, true, Si_pop, Ci_pop, Li_pop, X_pop, Y_pop, K, sim_data);
    % Save simulated data
    par_save(strcat('parsimdata_wsRPC_scen', num2str(scen),'_iter', num2str(sim_n)), sim_data);
    
    %% Sampling scenario 8: stratified random sample by subpop with equal allocation
    scen = 8;               
    n_s = [1000, 1000, 1000, 1000];  % Sample sizes for each subpop    
    % Obtain indices, weights, and data for sampled individuals
    sim_data = sample_indivs(N_s, n_s, S, true, Si_pop, Ci_pop, Li_pop, X_pop, Y_pop, K, sim_data);
    % Save simulated data
    par_save(strcat('parsimdata_wsRPC_scen', num2str(scen),'_iter', num2str(sim_n)), sim_data);
end




%% Testing code
% num_sims = 1;  
% N_s = [200, 200, 200, 200]; 
% n_s = 40;
% n_s = [10, 10, 10, 10];