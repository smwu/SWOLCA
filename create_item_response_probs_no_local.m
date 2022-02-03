% create_item_response_probs creates the true item response probabilities
% for each item-level combo, for all global and local classes
% Inputs:
%   p: Number of food items
%   d: Number of response levels for food items
%   clust_mode: Probability of true response level occurring
%   non_mode: Probability of other response levels occurring
%   sim_data: Structural array with the following fields:
%       true_pi: vector of true global class membership probabilities; Kx1
%       true_global_patterns: matrix of true item consumption levels for each global class; pxK
%   K: Number of global classes
% Outputs: returns updated sim_data with the following additional fields:
%   true_global_thetas: array of item response probs for global classes; pxKxd
function sim_data = create_item_response_probs_no_local(p, d, clust_mode, non_mode, sim_data, K)
    sim_data.true_global_thetas = ones(p, K, d) * non_mode;   % Initialize global item response probs to non_mode 

    for j = 1:p          % For each item
        for k = 1:K      % For each global class
            % Set prob true level occurs to clust_mode
            sim_data.true_global_thetas(j, k, sim_data.true_global_patterns(j, k)) = clust_mode; 
        end   
    end    
end 