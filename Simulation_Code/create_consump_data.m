% create_consump_data creates consumption data for each food item for all 
% individuals in the population.
% Inputs:
%   i_ind: Subject indicator
%   j_ind: Item indicator
%   s_ind: Subpopulation indicator
%   c_ind: Global class indicator
%   l_ind: Local class indicator
%   sim_data: Structural array with the following fields:
%       true_pi: vector of true global class membership probabilities; Kx1
%       true_lambda: vector of true local class membership probabilities; Lx1
%       true_global_patterns: matrix of true item consumption levels for each global class; pxK
%       true_local_patterns: array of true item consumption levels for each local class; pxSxL
%       true_global_thetas: array of item response probs for global classes; pxKxd
%       true_local_thetas: array of item response probs for local classes; pxSxLxd
%       nu: vector of subpop-specific probs of global assigment for all items; Sx1
%       Gsj: matrix of true assignments to global/local for each subpop and item; Sxp
%   X_pop: matrix of population consumption data for all indivs and items; Nxp
% Outputs:
%   X_pop: matrix of population consumption data for all indivs and items; Nxp
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