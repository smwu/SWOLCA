% create_item_response_probs creates the true item response probabilities
% for each item-level combo, for all global and local classes
% Inputs:
%   p: Number of food items
%   d: Number of response levels for food items
%   clust_mode: Probability of true response level occurring
%   non_mode: Probability of other response levels occurring
%   sim_data: Structural array with the following fields:
%       true_pi: vector of true global class membership probabilities; Kx1
%       true_lambda: vector of true local class membership probabilities; Lx1
%       true_global_patterns: matrix of true item consumption levels for each global class; pxK
%       true_local_patterns: array of true item consumption levels for each local class; pxSxL
%   S: Number of subpopulations
%   K: Number of global classes
%   L: Number of local classes per subpop
% Outputs: returns updated sim_data with the following additional fields:
%   true_global_thetas: array of item response probs for global classes; pxKxd
%   true_local_thetas: array of item response probs for local classes; pxSxLxd
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