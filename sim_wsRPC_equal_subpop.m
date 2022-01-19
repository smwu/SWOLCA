%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulated data for weighted supervised OFMM and RPC        %
% Programmer: SW                                             %
% Scenario 1: Equal subpop sizes                             %
%             Sample 100% from each subpop (full population) %
%             All weights equal to 1                         %     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data description:
% We assume individuals come from 4 subpopulations. Sampling of 4000 
% subjects is stratified by subpopulation. There are 3 global dietary 
% patterns, and for each subpopulation, there are 2 local dietary patterns. 
% Each food item for an individual belongs to a global pattern or to a 
% subpopulation-specific local pattern.

clear; clc;
rng(53188)     % Set seed

%% Set scenario specifications
num_sims = 100;                             % Number of simulation sets      
p = 50;                                     % Number of food items
d = 4;                                      % Number of response levels (assumed constant across items)
S = 4;                                      % Number of subpops
K = 3;                                      % Number of global classes
L = 2;                                      % Number of local classes per subpop
N_s = [20000, 20000, 20000, 20000];         % Subpopulation sizes           
N = sum(N_s);                               % Population size

sim_data.true_pi = [0.33 0.33 0.34];        % Global class membership proportions
sim_data.true_lambda = [0.5 0.5];           % Local class membership proportions
sim_data.true_xi = [1 1.1 1.2 0.9 1 -2 -1]; % True probit model coefficients
sim_data.nu = [0.5 0.4 0.7 0.8];            % Subpop-specific prob of global assigment for all items
clust_mode = 0.85;                          % Probability of true consumption level occurring
non_mode = 0.05;                            % Probability of other consumption levels occurring

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

%% Create and save simulations
for sim_n = 1:num_sims  % For each simulation iteration    
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
    
    %% Sampling scenario 1: full population
    scen = 1;  
    n_s = N;  % Sample size
    % Obtain indices, weights, and data for sampled individuals
    sim_data = sample_indivs(N, n_s, S, false, Si_pop, Ci_pop, Li_pop, X_pop, Y_pop, K, sim_data);
    % Save simulated data
    save(strcat('simdata_wsRPC_scen', num2str(scen),'_iter', num2str(sim_n)), 'sim_data');
    
    %% Sampling scenario 2: simple random sample
    scen = 2;                    
    n_s = 4000;  % Sample size 
    % Obtain indices, weights, and data for sampled individuals
    sim_data = sample_indivs(N, n_s, S, false, Si_pop, Ci_pop, Li_pop, X_pop, Y_pop, K, sim_data);
    % Save simulated data
    save(strcat('simdata_wsRPC_scen', num2str(scen),'_iter', num2str(sim_n)), 'sim_data');
    
    %% Sampling scenario 3: stratified random sample by subpop with proportional allocation
    scen = 3;                     
    n_s = round(0.05 .* N_s);  % Sample sizes for each subpop
    % Obtain indices, weights, and data for sampled individuals
    sim_data = sample_indivs(N_s, n_s, S, true, Si_pop, Ci_pop, Li_pop, X_pop, Y_pop, K, sim_data);
    % Save simulated data
    save(strcat('simdata_wsRPC_scen', num2str(scen),'_iter', num2str(sim_n)), 'sim_data');
    
    %% Sampling scenario 4: stratified random sample by subpop with equal allocation
    scen = 4;               
    n_s = [1000, 1000, 1000, 1000];  % Sample sizes for each subpop     
    % Obtain indices, weights, and data for sampled individuals
    sim_data = sample_indivs(N_s, n_s, S, true, Si_pop, Ci_pop, Li_pop, X_pop, Y_pop, K, sim_data);
    % Save simulated data
    save(strcat('simdata_wsRPC_scen', num2str(scen),'_iter', num2str(sim_n)), 'sim_data');
end




%% LOCAL FUNCTIONS

% create_item_response_probs creates the true item response probabilities
% for each item-level combo, for all global and local classes
% Inputs:
%   p: Number of food items
%   d: Number of response levels for food items
%   clust_mode: Probability of true response level occurring
%   non_mode: Probability of other response levels occurring
%   sim_data: Structural array with the following fields:
%       true_pi: Kx1 vector of true global class membership probabilities
%       true_lambda: Lx1 vector of true local class membership probabilities
%       true_global_patterns: pxK matrix of true item consumption levels for each global class
%       true_local_patterns: pxSxL array of true item consumption levels for each local class
%   S: Number of subpopulations
%   K: Number of global classes
%   L: Number of local classes per subpop
% Outputs: returns updated sim_data with the following additional fields:
%   true_global_thetas: pxKxd array of item response probs for global classes
%   true_local_thetas: pxSxLxd array of item response probs for local classes
function sim_data = create_item_response_probs(p, d, clust_mode, non_mode, sim_data, S, K, L)
    sim_data.true_global_thetas = ones(p, K, d) * non_mode;   % Initialize global item response probs to non_mode 
    sim_data.true_local_thetas = ones(p, S, L, d) * non_mode; % Initialize local item response probs to non_mode 

    for j = 1:p          % For each item
        for k = 1:K      % For each global class
            % Set prob true level occurs to clust_mode
            sim_data.true_global_thetas(j, k, sim_data.true_global_patterns(j, k)) = clust_mode; 
        end   
        for s = 1:S      % For each subpop
            for l = 1:L  % For each local class within the subpop
                % Set prob true level occurs to clust_mode 
                sim_data.true_local_thetas(j, s, l, sim_data.true_local_patterns(j, s, l)) = clust_mode; 
            end
        end
    end    
end 

% create_consump_data creates consumption data for each food item for all 
% individuals in the population.
% Inputs:
%   i_ind: Subject indicator
%   j_ind: Item indicator
%   s_ind: Subpopulation indicator
%   c_ind: Global class indicator
%   l_ind: Local class indicator
%   sim_data: Structural array with the following fields:
%       true_pi: Kx1 vector of true global class membership probabilities
%       true_lambda: Lx1 vector of true local class membership probabilities
%       true_global_patterns: pxK matrix of true item consumption levels for each global class
%       true_local_patterns: pxSxL array of true item consumption levels for each local class
%       true_global_thetas: pxKxd array of item response probs for global classes
%       true_local_thetas: pxSxLxd array of item response probs for local classes
%       nu: Sx1 vector of subpop-specific probs of global assigment for all items
%       Gsj: Sxp matrix of true assignments to global/local for each subpop and item
%   X_pop: Nxp matrix of population consumption data for all indivs and items
% Outputs:
%   X_pop: Nxp matrix of population consumption data for all indivs and items
function X_pop = create_consump_data(i_ind, j_ind, s_ind, c_ind, l_ind, sim_data, X_pop)

    if sim_data.Gsj(s_ind, j_ind) == 1  % Item assigned to global pattern
        % Draw from a Multinomial to obtain consump level based on global item response probs 
        item_resp = sim_data.true_global_thetas(j_ind, c_ind, :);
        draw = mnrnd(1, item_resp(:));
        X_pop(i_ind, j_ind) = find(draw == 1);

    else  % Item assigned to local pattern
        % Draw from a Multinomial to obtain consump level based on local item response probs
        item_resp = sim_data.true_local_thetas(j_ind, s_ind, l_ind, :);
        draw = mnrnd(1, item_resp(:)); 
        X_pop(i_ind, j_ind) = find(draw == 1);
    end
    
end

% sample_indivs takes in subpop specs and outputs indices and weights for 
% sampled individuals.
% Inputs:
%   N_s: If stratified, Sx1 vector of subpop sizes; else, population size
%   n_s: If stratified, Sx1 vector of subpop sample sizes; else, sample size
%   S: Number of subpops
%   strat: Boolean specifying whether to stratify sampling by subpop
%   Si_pop: Nx1 vector of subpop assignments for all indivs
%   Ci_pop: Nx1 vector of global class assignments for all indivs
%   Li_pop: Nx1 vector of local class assignments for all indivs
%   X_pop: Nxp matrix of population consumption data for all indivs and items
%   Y_pop: Nx1 vector of population outcome data for all indivs
%   K: Number of global classes
%   sim_data: Structural array with the following fields:
%       true_pi: Kx1 vector of true global class membership probabilities
%       true_lambda: Lx1 vector of true local class membership probabilities
%       true_global_patterns: pxK matrix of true item consumption levels for each global class
%       true_local_patterns: pxSxL array of true item consumption levels for each local class
%       true_global_thetas: pxKxd array of item response probs for global classes
%       true_local_thetas: pxSxLxd array of item response probs for local classes
%       nu: Sx1 vector of subpop-specific probs of global assigment for all items
%       Gsj: Sxp matrix of true assignments to global/local for each subpop and item
%       true_xi: (K+S)x1 vector of true probit model coefficients
%       true_Phi: Nx1 vector of true probit means, P(y_i=1|Q)
% Outputs: Updated sim_data with the following additional fields:
%   samp_ind: nx1 vector of indices of sampled indivs
%   sample_wt: nx1 vector of sampling weights for sampled indivs
%   norm_const: Sx1 vector of normalization constants, one per subpop  
%   true_Si: nx1 vector of subpop assigns for sampled indivs
%   true_Ci: nx1 vector of global class assigns for sampled indivs
%   true_Li: nx1 vector of local class assigns for sampled indivs
%   X_data: nxp matrix of consumption data for sampled indivs and items
%   Y_data: nx1 vector of outcome data for sampled indivs and items
%   true_K: Number of global classes
function sim_data = sample_indivs(N_s, n_s, S, strat, Si_pop, Ci_pop, Li_pop, X_pop, Y_pop, K, sim_data) 
    if strat == false                                        % If not stratifying by subpop
        sim_data.samp_ind = randsample(N_s, n_s);            % Indices of randomly sampled indivs
        wt_s = N_s / n_s;                                    % Sampling weight
        pop_wt = ones(N_s, 1) * wt_s;                        % Weights, temp applied to all popn indivs
        sim_data.sample_wt = pop_wt(sim_data.samp_ind);      % Weights for sampled indivs 
        sim_data.norm_const = sum(sim_data.sample_wt / N_s); % Normalization constant. Default is 1

    else 
        sim_data.samp_ind = cell(S, 1);   % Initialize indices of sampled indivs, grouped by subpop
        sim_data.sample_wt = cell(S, 1);  % Initialize sampling weights for sampled indivs, grouped by subpop
        sim_data.norm_const = cell(S, 1); % Initialize normalization constant for each subpop
      
        % Obtain indices and weights for sampled individuals in each subpop
        for s = 1:S  % For each subpopulation
            % Indices of randomly sampled indivs from subpop s
            sim_data.samp_ind{s} = randsample(N_s(s), n_s(s));     
            wt_s = N_s(s) / n_s(s);             % Sampling weight for subpop s
            pop_wt_s = ones(N_s(s), 1) * wt_s;  % Weights for subpop s, temp applied to all popn indivs
            % Weights for sampled indivs in subpop s
            sim_data.sample_wt{s} = pop_wt_s(sim_data.samp_ind{s});        
            % Normalization constant for subpop s weights
            sim_data.norm_const{s} = sum(sim_data.sample_wt{s}) / N_s(s); 
        end
        
        % Convert cell arrays to vectors through vertical concat 
        sim_data.samp_ind = vertcat(sim_data.samp_ind{:});  
        sim_data.sample_wt = vertcat(sim_data.sample_wt{:});
        sim_data.norm_const = vertcat(sim_data.norm_const{:});
    end    

    % Obtain observed sample data
    sim_data.true_Si = Si_pop(sim_data.samp_ind);   % Subpop assignments for sampled indivs
    sim_data.true_Ci = Ci_pop(sim_data.samp_ind);   % True global class assignments for sampled indivs
    sim_data.true_Li = Li_pop(sim_data.samp_ind);   % True local class assignments for sampled indivs
    sim_data.X_data = X_pop(sim_data.samp_ind, :);  % Consumption data for sampled indivs
    sim_data.Y_data = Y_pop(sim_data.samp_ind);     % Outcome data for sampled indivs
    sim_data.true_K = K;                            % Rename variable for number of classes
end

%% Testing code
% num_sims = 1;  
% N_s = [200, 200, 200, 200]; 
% n_s = 40;
% n_s = [10, 10, 10, 10];