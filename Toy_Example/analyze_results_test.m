% analyze_results obtains reduced number of classes, posterior estimates of 
% the params, and analysis metrics
% Inputs: 
%   MCMC_out: output from run_MCMC_latent or wtd_run_MCMC_latent function
%   post_MCMC_out: output from post_process function
%   data_vars: output from wtd_get_data_vars or get_data_vars function
%   q_dem: Matrix of demographic covariates to include in regression, in cell-means format
%   S: number of demographic covariates in the probit model
%   p_cov: Number of covariates in probit model
% Outputs: returns analysis structural array with the following fields:
%   k_red: Reduced number of classes after removing duplicates
%   pi_red: Matrix of posterior class membership probs for all MCMC runs; Mx(k_red)
%   theta_red: Array of posterior iter-response probs for all MCMC runs; Mxpx(k_red)p(d_max)
%   xi_red: Matrix of posterior probit model coefs for all MCMC runs; Mx(S+k_red)
%   pi_med: Vector of posterior median estimates across MCMC runs for class membership probs; 1x(k_red)
%   theta_med: Array of posterior median estimates across MCMC runs for item-response probs; px(k_red)x(d_max)
%   xi_med: Vector of posterior median estimates across MCMC runs for probit model coefs; 1x(S+k_red)
%   z_i: Vector of latent probit variable for each individual for the last MCMC run; nx1
%   max_prob_ci: Vector of max posterior class assignment probs for each individual; nx1
%   c_i: Vector of class assignment for each individual using posterior median estimates; nx1
%   Phi_med: estimated probit model mean using posterior median estiates; nx1 
%   loglik_med: Log-lik using posterior median estimates of all params; nx1
%   dic6: Assesses model goodness-of-fit with adaptation that penalizes overfitting
%   aebic: Assesses model goodness-of-fit with asymp consistent and reparam-invariant metric 
function analysis = analyze_results_test(MCMC_out, post_MCMC_out, data_vars, q_dem, S, p_cov)
    % Identify unique classes     
    % Array of posterior median item-response probs
    theta_med_temp = reshape(median(post_MCMC_out.theta), [data_vars.p, post_MCMC_out.k_med, data_vars.d_max]);  
    % Matrix of most likely consump level based on item-response probs, for each item and class across MCMC runs
    [~, ind0] = max(theta_med_temp, [], 3); 
    t_ind0 = transpose(ind0);
    [~, ia, ~] = unique(t_ind0, 'rows');  % Vector of class indices corresponding to unique classes
    analysis.k_red = length(ia);          % Number of unique classes
   
    % Posterior MCMC parameter chains for the unique classes, re-normalized
    analysis.pi_red = post_MCMC_out.pi(:, ia);
    % analysis.pi_red = analysis.pi_red ./ sum(analysis.pi_red, 2);
    analysis.theta_red = post_MCMC_out.theta(:, :, ia, :);
    % analysis.theta_red = analysis.theta_red ./ sum(analysis.theta_red, 4);
    analysis.xi_red = post_MCMC_out.xi(:, [1:S (S+ia')]);
    
    % Posterior median estimates for the unique classes
    analysis.pi_med = median(analysis.pi_red);  % Class probs 
    analysis.pi_med = analysis.pi_med / sum(analysis.pi_med);
    % Item-response probs 
    analysis.theta_med = reshape(median(analysis.theta_red), [data_vars.p, analysis.k_red, data_vars.d_max]);  
    analysis.theta_med = analysis.theta_med ./ sum(analysis.theta_med, 3);
    analysis.xi_med = median(analysis.xi_red);  % Probit coefs   
    
    analysis.z_i = MCMC_out.z_i(end, :)'; % Obtain last values of probit latent variable

    % Update class membership probit component using unique classes and updated posterior median estimates for xi
    indiv_lik_probit_class = zeros(data_vars.n, analysis.k_red); % Initialize indiv probit lik matrix for each indiv and class
    for k = 1:analysis.k_red
       q_class = zeros(data_vars.n, analysis.k_red);            % Initialize design matrix for class membership
       q_class(:, k) = ones(data_vars.n, 1);                    % Temporarily assign all indivs to class k
       Q_temp = [q_dem q_class];                                % Form temporary covariate design matrix
       Phi_temp = normcdf(Q_temp * transpose(analysis.xi_med)); % Vector of P(y_i=1|Q)
       Phi_temp(Phi_temp == 1) = 1 - 1e-10;                     % Adjust extremes
       Phi_temp(Phi_temp == 0) = 1e-10;
       % Update indiv probit log-lik for class k. Helps avoid rounding issues
       loglik_probit = data_vars.y .* log(Phi_temp) + (1-data_vars.y) .* log(1-Phi_temp);
       indiv_lik_probit_class(:, k) = exp(loglik_probit);    % Update indiv probit likelihood for class k
    end  
    
    % Update class membership, c_i, using unique classes and updated posterior median estimates for pi and theta
 	prob_ci_numer = zeros(data_vars.n, analysis.k_red);  % Initialize indiv OFMM lik matrix for each indiv and class
    for k = 1:analysis.k_red
        % Matrix of unique item-level probs for class k
        theta_k = reshape(analysis.theta_med(:, k, :), data_vars.p, data_vars.d_max); 
        % Matrix of categ likelihood \prod_{r=1}^d \theta_{jr|k}^{I(x_ij=r)} for each indiv and item 
        categ_lik = reshape(theta_k(data_vars.lin_idx), [data_vars.n, data_vars.p]); 
        % Update numerator for posterior P(c_i=k|-) (indiv likelihood)
        prob_ci_numer(:, k) = analysis.pi_med(k) * prod(categ_lik, 2) .* indiv_lik_probit_class(:, k);   
    end
    % Matrix of updated posterior P(c_i=k|-) for each indiv and class
    prob_ci = bsxfun(@times, prob_ci_numer, 1 ./ (sum(prob_ci_numer, 2)));
    % Max posterior class assign probs and class assigns for each indiv using posterior median estimates
    [analysis.max_prob_ci, analysis.c_i] = max(prob_ci, [], 2);  
%     x_ci = mnrnd(1, prob_ci);                   % Matrix of updated c_i drawn from Mult(1, prob_ci). Each row is draw for an indiv
%     [row, col] = find(x_ci);                    % Row and col indices of each nonzero element of x_ci; col is class assignment for each indiv
%     sorted = sortrows([row col], 1);            % Sort indices in ascending row index order, to match subject index
%     analysis.c_i = sorted(:, 2);                % Vector of updated class assignment, c_i, for each individual
      
    % Obtain probit model mean using median parameter estimates
%    Q_med = [q_dem x_ci];                             % Design matrix with updated class memberships using median estimates   
    Q_med = [q_dem dummyvar(analysis.c_i)];  % Design matrix with updated class memberships using median estimates 
    analysis.Phi_med = normcdf(Q_med * transpose(analysis.xi_med)); % Probit model mean using median estimates
    analysis.Phi_med(analysis.Phi_med == 0) = 1e-10;  % Adjust extremes
    analysis.Phi_med(analysis.Phi_med == 1) = 1 - 1e-10;
    
    % Indiv observed log likelihood calculated using posterior median param estimates
    analysis.loglik_med = log(sum(prob_ci_numer, 2)); 
    
    % DIC is a metric for model goodness of fit. It consists of two terms:
        % 1) median(loglik): calculate log-lik at each MCMC iter then get the median
        % 2) sum(loglik_med): get posterior estimates of each param, then use these to calculate the log-lik
    % DIC-6 is an adaptation that penalizes overfitting
    analysis.dic6 = -6 * median(post_MCMC_out.loglik) + 4 * sum(analysis.loglik_med); 
    
    % AEBIC is a metric for model goodness of fit that is asymptotically
    % consistent and invariant to reparameterization
    num_params = numel(analysis.pi_med) + numel(analysis.theta_med) + numel(analysis.xi_med);
    gamma = 1;
    t1 = -2 * data_vars.n * mean(post_MCMC_out.loglik ./ data_vars.n); % -2*n*E[1/n*sum(indiv_log_lik)]
    t2 = num_params * log(data_vars.n);                                % |S|*log(n)
    t3 = 2 * gamma * num_params * log(p_cov);                          % 2*gamma*|S|*log(p)
    analysis.aebic = t1 + t2 + t3;
end