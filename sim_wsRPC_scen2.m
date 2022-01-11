%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulated data for weighted supervised OFMM and RPC        %
% Programmer: SW                                             %
% Scenario 2: Unequal subpop sizes                           %
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
num_sims = 1;  % Number of simulation sets
scen = 2;      % Scenario
rng(53188)     % Set seed

for sim_n = 1:num_sims                     % For each simulation iteration
    %% Set scenario specification
    S = 4;                              % Number of subpops
    K = 3;                              % Number of global classes
    L = 2;                              % Number of local classes per subpop
    N_s = [500, 1250, 1375, 875];       % Subpop sizes
    n_s = N_s;                          % Sampled subpop sizes
    n = sum(n_s);                       % Sample size
    p = 50;                             % Number of food items
    d = 4;                              % Number of response levels (assumed constant across items)
    
    %% Create subpop assignments, true global and local class assignments, and sampling weights for each subject
    subpop = cell(S, 1);      % Initialize subpop assignments for all indivs
    subpop_samp = cell(S, 1); % Initialize subpop assignments for sampled indivs    
    true_Ci = cell(S, 1);     % Initialize true global class assignments for sampled indivs
    true_Li = cell(S, 1);     % Initialize true local class assignments for sampled indivs    
    sample_wt = cell(S, 1);   % Initialize sampling weights for sampled indivs, grouped by subpop
    norm_const = cell(S, 1);  % Initialize normalization constant for each subpop. Default is 1
    
    for s = 1:S                                % For each subpopulation
        subpop{s} = ones(N_s(s), 1) * s;       % Assign subpop for popn indivs
        subpop_samp{s} = ones(n_s(s), 1) * s;  % Assign subpop for sampled indivs        
        % Create global and local class assigns and weights for sampled indivs in subpop s
        [true_Ci{s}, true_Li{s}, sample_wt{s}, norm_const{s}] = create_classes_weights(N_s(s), n_s(s), K, L);
    end
    
    subpop = vertcat(subpop{:});  % Convert cell array to vector through vertical concat
    subpop_samp = vertcat(subpop_samp{:});
    true_Ci = vertcat(true_Ci{:});  
    true_Li = vertcat(true_Li{:});
    sample_wt = vertcat(sample_wt{:});
    norm_const = vertcat(norm_const{:});
    
    %% Create true comsumption patterns for each diet profile class      
    global1 = [ones(0.5 * p, 1) * 3;  % Global profile patterns
               ones(0.5 * p, 1) * 1];
    global2 = [ones(0.2 * p, 1) * 2; 
               ones(0.8 * p, 1) * 4];
    global3 = [ones(0.2 * p, 1) * 1; 
               ones(0.4 * p, 1) * 2; 
               ones(0.4 * p, 1) * 3];
    % p x K matrix of true consumption levels for each global class
    true_global_patterns = [global1 global2 global3];  
       
    local11 = ones(p, 1) * 1;  % Local profile patterns
    local12 = [repmat([1; 2; 3; 4], [12, 1]); 1; 1];
    local21 = ones(p, 1) * 2;
    local22 = [repmat([2; 4; 1; 3], [12, 1]); 2; 2];
    local31 = ones(p, 1) * 3;
    local32 = [repmat([3; 1; 4; 2], [12, 1]); 3; 3];
    local41 = ones(p, 1) * 4;
    local42 = [repmat([4; 3; 2; 1], [12, 1]); 4; 4];
    % p x S x L array of true consumption levels for each local class
    true_local_patterns = [local11 local21 local31 local41];
    true_local_patterns(:,:,2) = [local12 local22 local32 local42];
        
    %% Create true item response probabilities
    clust_mode = 0.85;  % Probability of true level occurring
    non_mode = 0.05;    % Probability of other levels occurring
    [global_thetas, local_thetas, global_thetas_sum, local_thetas_sum] = create_item_response_probs(p, d, clust_mode, non_mode, true_global_patterns, true_local_patterns, S, K, L);
    
    %% Create deviation indicator for global or local assignment
    % Assume all indivs in the same subpop have the same items that are global/local, i.e., G_ij = G_sj
    % Assume all items for those in the same subpop have the same deviation prob, i.e., nu_sj = nu_s
    
    nu = [0.5 0.4 0.7 0.8];  % Subpop-specific prob of global assigment for all items
    true_G = zeros(S, p);    % Initialize matrix of true global/local assign for each subpop and item
    for s = 1:S              % For each subpop
        % Randomly sample vector of global/local assigns from Ber(nu(s))
        true_G(s, :) = binornd(1, nu(s), [1, p]);  
    end
    nu_sum = sum(true_G, 2);  % Number of global items for each subpop
    
    %% Create observed consumption data
    rand_unif = rand(n, p);     % Matrix of Unif(0,1) random variables
    sample_data = zeros(n, p);  % Initialize sample observed data for each individual and item
    
    for i = 1:n
        for j = 1:p
            % Create observed consumption data for each indiv and item
            sample_data = create_consump_data(i, j, subpop_samp(i), true_Ci(i), true_Li(i), d, true_G, rand_unif, global_thetas_sum, local_thetas_sum, sample_data);
        end
    end  
    
    %% Create true probit model and parameters        
    true_xi = [1 1.1 1.2 0.9 1 -2 -1];  % True probit model coefficients
    q_dem = zeros(n, S);                % Initialize design matrix of dem covs (subpop) in cell-means format
    for s = 1:S                         % For each subpop
        q_dem(subpop_samp == s, s) = 1; % For all indivs in subpop s, set s-th column to 1
    end
    q_class = zeros(n, K);              % Initialize design matrix for global class membership in cell-means format
    for k = 1:K                         % For each global class
        q_class(true_Ci == k, k) = 1;   % For all indivs in global class k, set k-th column to 1
    end
    
    Q = [q_dem q_class];                    % Design matrix with dem covs and global classes
    true_lin_pred = Q * transpose(true_xi); % True linear predictor for sampled indivs. Mean of truncated normal dist
    true_Phi = normcdf(true_lin_pred);      % True probit mean, P(y_i=1|Q)
    true_y = binornd(1, true_Phi);          % True outcome for sampled indivs
    
    % Save simulated data
    true_K = K;
    save(strcat('simdata_wsRPC_scen', num2str(scen),'_iter', num2str(sim_n)), 'subpop_samp', 'true_Ci', 'true_Li', 'sample_wt', 'norm_const', 'true_global_patterns', 'true_local_patterns', 'global_thetas', 'local_thetas', 'nu', 'true_G', 'sample_data', 'true_xi', 'true_Phi', 'true_y', 'true_K');
end




%% LOCAL FUNCTIONS

% create_classes_weights takes in subpop specs and outputs global and local
% class assignments, as well as sampling weights, for each sampled subject
% Inputs:
%   N_s_ind: Population size of subpop s
%   n_s_ind: Sample size of subpop s
%   K: Number of global classes
%   L: Number of local classes per subpop
% Outputs:
%   true_Ci_s: vector of true global class assigns for sampled indivs in subpop s
%   true_Li_s: vector of true local class assigns for sampled indivs in subpop s
%   sample_wt_s: vector of sampling weights for sampled indivs in subpop s
%   norm_const_s: normalization constant for subpop s
function [true_Ci_s, true_Li_s, sample_wt_s, norm_const_s] = create_classes_weights(N_s_ind, n_s_ind, K, L) 
    C_i = randi([1 K], N_s_ind, 1);  % Randomly generate global class assigns for subpop s 
    L_i = randi([1 L], N_s_ind, 1);  % Randomly generate local class assigns for subpop s 
    
    samp_ind = randsample(N_s_ind, n_s_ind); % Indices of randomly sampled indivs from subpop s
    true_Ci_s = C_i(samp_ind);               % True global class assigns for sampled indivs in subpop s
    true_Li_s = L_i(samp_ind);               % True local class assigns for sampled indivs in subpop s
    
    wt_s = N_s_ind / n_s_ind;                  % Sampling weight for subpop s
    pop_wt_s = ones(N_s_ind, 1) * wt_s;        % Weights for subpop s, temp applied to all popn indivs
    sample_wt_s = pop_wt_s(samp_ind);          % Weights for sampled indivs in subpop s
    norm_const_s = sum(sample_wt_s) / N_s_ind; % Normalization constant for subpop s weights
end

% create_item_response_probs creates the true item response probabilities
% for each item-level combo, for all global and local classes
% Inputs:
%   p: Number of food items
%   d: Number of response levels for food items
%   clust_mode: Probability of true response level occurring
%   non_mode: Probability of other response levels occurring
%   true_global_patterns: true item consumption levels for each global class
%   true_local_patterns: true item consumption levels for each local class
%   S: Number of subpopulations
%   K: Number of global classes
%   L: Number of local classes per subpop
% Outputs:
%   global_thetas: Item response probs for global classes
%   local_thetas: Item response probs for local classes
%   global_thetas_sum: Cumulative sums of item response probs for global classes
%   local_thetas_sum: Cumulative sums of item response probs for local classes
function [global_thetas, local_thetas, global_thetas_sum, local_thetas_sum] = create_item_response_probs(p, d, clust_mode, non_mode, true_global_patterns, true_local_patterns, S, K, L)
    global_thetas = ones(p, d, K) * non_mode;   % Initialize global item response probs to non_mode 
    local_thetas = ones(p, d, S, L) * non_mode; % Initialize local item response probs to non_mode 
    global_thetas_sum = zeros(p, d + 1, K);     % Initialize cumsum of global item response probs
    local_thetas_sum = zeros(p, d + 1, S, L);   % Initialize cumsum of local item response probs

    for j = 1:p          % For each item
        for k = 1:K      % For each global class
            % Set prob true level occurs to clust_mode
            global_thetas(j, true_global_patterns(j, k), k) = clust_mode; 
            % Set cumsum of global item response probs for item j in class k
            global_thetas_sum(j, :, k) = [0 cumsum(global_thetas(j, :, k), 2)];
        end   
        for s = 1:S      % For each subpop
            for l = 1:L  % For each local class within the subpop
                % Set prob true level occurs to clust_mode 
                local_thetas(j, true_local_patterns(j, s, l), s, l) = clust_mode; 
                % Set cumsum of local item response probs for item j in subpop s class l
                local_thetas_sum(j, :, s, l) = [0 cumsum(local_thetas(j, :, s, l), 2)];
            end
        end
    end    
end 

% create_consump_data creates observed consumption data for each food item for all sampled individuals
% Inputs:
%   i_ind: Subject indicator
%   j_ind: Item indicator
%   s_ind: Subpopulation indicator
%   k_ind: Global class indicator
%   l_ind: Local class indicator
%   d: Number of item response levels
%   true_G: Matrix of true assignments to global/local for each subpop and item
%   rand_unif: Matrix of Unif(0,1) random variables
%   global_thetas_sum: Cumulative sums of item response probs for global classes
%   local_thetas_sum: Cumulative sums of item response probs for local classes
%   sample_data: Matrix of sample observed data for each individual and item
% Outputs:
%   sample_data: Updated matrix of sample observed data for each individual and item
function sample_data = create_consump_data(i_ind, j_ind, s_ind, k_ind, l_ind, d, true_G, rand_unif, global_thetas_sum, local_thetas_sum, sample_data)
    if true_G(s_ind, j_ind) == 1  % Item assigned to global pattern
        % Randomly assign a consumption level based on item response probs
        % Equiv to drawing from a Categorical dist with probs given by the thetas
        for r = 1:d
            ind = rand_unif(i_ind, j_ind) > global_thetas_sum(j_ind, r, k_ind) & rand_unif(i_ind, j_ind) <= global_thetas_sum(j_ind, r+1, k_ind);
            if ind == 1
                sample_data(i_ind, j_ind) = r;
            end   
        end
    else  % Item assigned to local pattern
        for r = 1:d
            ind = rand_unif(i_ind, j_ind) > local_thetas_sum(j_ind, r, s_ind, l_ind) & rand_unif(i_ind, j_ind) <= local_thetas_sum(j_ind, r+1, s_ind, l_ind);
            if ind == 1
                sample_data(i_ind, j_ind) = r;
            end   
        end
    end
end
