% analyze_results_unsup_latent obtains reduced number of classes, posterior 
% estimates of the params, and analysis metrics for the wOFMM model
% Inputs:
%   post_MCMC_out: output from post_process_latent function
%   data_vars: output from wtd_get_data_vars_latent or get_data_vars_latent function
% Outputs: returns analysis structural array with the following fields:
%   k_red: Reduced number of classes after removing duplicates
%   pi_red: Matrix of posterior class membership probs for all MCMC runs; Mx(k_red)
%   theta_red: Array of posterior iter-response probs for all MCMC runs; Mxpx(k_red)p(d_max)
%   pi_med: Vector of posterior median estimates across MCMC runs for class membership probs; 1x(k_red)
%   theta_med: Array of posterior median estimates across MCMC runs for item-response probs; px(k_red)x(d_max)
%   max_prob_ci: Vector of max posterior class assignment probs for each individual; nx1
%   c_i: Vector of class assignment for each individual using posterior median estimates; nx1
function analysis = analyze_results_unsup_latent(post_MCMC_out, data_vars)
    % Identify unique classes   
    % Array of posterior median item-response probs
    theta_med_temp = reshape(median(post_MCMC_out.theta), [data_vars.p, post_MCMC_out.k_med, data_vars.d_max]);  
    % Matrix of most likely consump level based on item-response probs, for each item and class across MCMC runs
    [~, ind0] = max(theta_med_temp, [], 3); 
    t_ind0 = transpose(ind0);
    [u_mat, ia, ~] = unique(t_ind0, 'rows');  % Vector of class indices corresponding to unique classes
    analysis.k_red = length(ia);          % Number of unique classes
        
    % Posterior MCMC chains for pi after combining duplicated classes, re-normalized
    analysis.pi_red = post_MCMC_out.pi(:, ia);    % Restrict to unique classes
    has_dupes = size(u_mat, 1) < size(t_ind0, 1); % Check for duplicates
    if has_dupes
        ind_dup = setdiff(1:size(t_ind0,1), ia);  % Find indices of duplicated classes
        for i = 1:length(ind_dup)                    % Iterate over dupe rows
            for j = 1:length(ia)                     % Iterate over non-dupe rows
                i_dup = ind_dup(i);               
                i_orig = ia(j);
                if t_ind0(i_dup, :) == t_ind0(i_orig, :)
                    % Combine pi for duplicated patterns
                    analysis.pi_red(:, i_orig) = analysis.pi_red(:, i_orig) + post_MCMC_out.pi(:, i_dup);
                    %disp(median(analysis.pi_red));
                end
            end
        end
    end  
    analysis.pi_red = analysis.pi_red ./ sum(analysis.pi_red, 2);
    
    % Posterior MCMC parameter chains for the unique classes, re-normalized
    analysis.theta_red = post_MCMC_out.theta(:, :, ia, :);
    analysis.theta_red = analysis.theta_red ./ sum(analysis.theta_red, 4);   

    % Posterior median estimates for the unique classes
    analysis.pi_med = median(analysis.pi_red);  % Class probs 
    analysis.pi_med = analysis.pi_med / sum(analysis.pi_med);
    % Item-response probs 
    analysis.theta_med = reshape(median(analysis.theta_red), [data_vars.p, analysis.k_red, data_vars.d_max]);  
    analysis.theta_med = analysis.theta_med ./ sum(analysis.theta_med, 3);    
    
    % Update class membership, c_i, using unique classes and updated posterior median estimates for pi and theta
 	prob_ci_numer = zeros(data_vars.n, analysis.k_red);  % Initialize indiv OFMM lik matrix for each indiv and class
    for k = 1:analysis.k_red
        % Matrix of unique item-level probs for class k
        theta_k = reshape(analysis.theta_med(:, k, :), data_vars.p, data_vars.d_max); 
        % Matrix of categ likelihood \prod_{r=1}^d \theta_{jr|k}^{I(x_ij=r)} for each indiv and item 
        categ_lik = reshape(theta_k(data_vars.lin_idx), [data_vars.n, data_vars.p]); 
        % Update numerator for posterior P(c_i=k|-) (indiv OFMM likelihood)
        prob_ci_numer(:, k) = analysis.pi_med(k) * prod(categ_lik, 2);   
    end
    % Matrix of updated posterior P(c_i=k|-) for each indiv and class
    prob_ci = bsxfun(@times, prob_ci_numer, 1 ./ (sum(prob_ci_numer, 2)));
    % Max posterior class assign probs and class assigns for each indiv using posterior median estimates
    [analysis.max_prob_ci, analysis.c_i] = max(prob_ci, [], 2);   
    
%     %% Obtain probit model mean using median parameter estimates
%     % Update latent probit variable, z_i
%     % Factor variable coding design matrix with demographic covariates and initial latent classes
%     Q = zeros(data_vars.n, S * k_max);
%     for s = 1:S
%         for k = 1:k_max
%             Q(:, ((s-1) * k_max) + k) = q_dem(:, s) .* OFMM_params.x_ci(:, k);
%         end
%     end    
%     lin_pred = Q * transpose(probit_params.xi);  % Linear predictor, Q*xi. Mean of truncated normal dist
% 
% 
%     % Update posterior probit model coefficients, xi, incorporating normalized weights                 
%     diag_ind = 1:data_vars.n;                                     % Indices to indicate matrix diagonal
%     W_tilde = sparse(diag_ind, diag_ind, data_vars.wt_kappa);     % Create sparse diagonal normalized weight matrix
%     Q_sparse = sparse(Q);                                         % Convert to sparse design matrix to save memory
%     Sig_post = inv(Sig_0) + full(Q_sparse' * W_tilde * Q_sparse); % Precision of posterior Normal dist for xi
%     mu_right = (Sig_0 \ mu_0) + (Q_sparse' * W_tilde * probit_params.z_i); % Right part of mean of posterior Normal dist for xi
%     mu_post = Sig_post \ mu_right;                                % Mean of posterior Normal dist for xi
%     probit_params.xi = mvnrnd(mu_post, inv(Sig_post));            % Draw from posterior dist: N(mu_post, Sig_post^{-1})  
%     
% 
%     
%     % Factor variable coding design matrix with updated class memberships using median estimates 
%     Q_med = zeros(data_vars.n, S * analysis.k_red);
%     x_ci = dummyvar(analysis.c_i);
%     for s = 1:S
%         for k = 1:analysis.k_red
%             Q_med(:, ((s-1) * analysis.k_red) + k) = q_dem(:, s) .* x_ci(:, k);
%         end
%     end
%     analysis.Phi_med = normcdf(Q_med * transpose(analysis.xi_med)); % Linear predictor, Q*xi. Probit model mean using median estimates
%     analysis.Phi_med(analysis.Phi_med == 0) = 1e-10;  % Adjust extremes
%     analysis.Phi_med(analysis.Phi_med == 1) = 1 - 1e-10;
    
    
end