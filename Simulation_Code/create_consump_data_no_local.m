% create_consump_data creates consumption data for each food item for all 
% individuals in the population.
% Inputs:
%   i_ind: Subject indicator
%   j_ind: Item indicator
%   c_ind: Global class indicator
%   sim_data: Structural array with the following fields:
%       true_pi: vector of true global class membership probabilities; Kx1
%       true_lambda: vector of true local class membership probabilities; Lx1
%       true_global_patterns: matrix of true item consumption levels for each global class; pxK
%       true_global_thetas: array of item response probs for global classes; pxKxd
%   X_pop: matrix of population consumption data for all indivs and items; Nxp
% Outputs:
%   X_pop: matrix of population consumption data for all indivs and items; Nxp
function X_pop = create_consump_data_no_local(i_ind, j_ind, c_ind, sim_data, X_pop)

    % Draw from a Multinomial to obtain consump level based on global item response probs 
    item_resp = sim_data.true_global_thetas(j_ind, c_ind, :);
    draw = mnrnd(1, item_resp(:));
    X_pop(i_ind, j_ind) = find(draw == 1);
    
end