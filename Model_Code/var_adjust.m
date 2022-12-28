function [analysis] = var_adjust(analysis, samp_data, data_vars, psu_indic, R)
    
    % Get posterior means
    pi_est = mean(analysis.pi_red);
    theta_est  = mean(analysis.theta_red);
    xi_est = mean(analysis.xi_red);

    % Estimate Hessian
    log_post = 
    H_hat = hessian(log_post, pi, theta, xi);
    
    j_r = zeros(R, 1);
    % For each replicate
    for r = 1:R
        x_r = cell(S, 1);
        y_r = cell(S, 1);
        wt_r = cell(S, 1);
        for s = 1:S  % for each stratum
            % sample half of the PSUs within stratum
            subset = find(samp_data.true_Si == s);                   % Subset to the indices of indivs in subpop s
            n_s = length(subset);
            ind_s = sort(randsample(subset, round(n_s/2))); % Randomly sample n_s indices from subpop s; sort by ascending index
            x_r{s} = samp_data.X_data(ind_s, :); 
            y_r{s} = samp_data.Y_data(ind_s);
            wt_s = data_vars.wt_kappa(ind_s);                        
            wt_r{s} = wt_s * 2;                                         % New sampling weight for subpop s 
        end
        % Combine samples across strata
        x_r = sort(vertcat(x_r{:}));
        y_r = sort(vertcat(y_r{:}));
        wt_r = sort(vertcat(wt_r{:}));

        % Normalize weights
        wt_r_kappa = wt_r * length(wt_r) / sum(wt_r);

        % Compute J_hat component by taking the gradient of the weighted log-posterior
        new_log_post = 
        j_r(r) = gradient(new_log_post)
    end

    J_hat = vcov(j_r);

    % Compute adjustment
    H_inv = inv(H_hat);
    V_1 = (H_hat \ J_hat) / H_hat;  % H_hat^{-1} * J_hat * H_hat^{-1}
    R_1 = chol(V_1);
    R_2_inv = chol(H_inv);
    R_2R_1 = R_2_inv \ R_1;  % R_2 * R_1

    % Adjust samples
    analysis.pi_adj = (analysis.pi_red - pi_est) * R_2R_1 + pi_est;
    analysis.theta_adj = (analysis.theta_red - theta_est) * R_2R_1 + theta_est;
    analysis.xi_adj = (analysis.xi_red - xi_est) * R_2R_1 + xi_est; 


end
