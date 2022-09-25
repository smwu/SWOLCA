% analyze_results_latent obtains reduced number of classes, posterior estimates of 
% the params, and analysis metrics with the latent variable formulation
% Inputs:
%   MCMC_out: output from run_MCMC_latent or wtd_run_MCMC_latent function
%   post_MCMC_out: output from post_process_latent function
%   data_vars: output from wtd_get_data_vars_latent or get_data_vars_latent function
%   q_dem: matrix of demographic covariates to include in regression, in cell-means format
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
function analysis = analyze_results_latent(MCMC_out, post_MCMC_out, data_vars, q_dem, S, p_cov)
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
    % Factor variable version
    analysis.xi_red = zeros(size(analysis.pi_red, 1), S * analysis.k_red);
    for s = 1:S
        analysis.xi_red(:, ((s-1)*analysis.k_red) + (1:analysis.k_red)) = post_MCMC_out.xi(:, (s-1)*post_MCMC_out.k_med + ia);
    end
    %         analysis.xi_red = post_MCMC_out.xi(:, [1:S (S+ia')]);     

    % Posterior median estimates for the unique classes
    analysis.pi_med = median(analysis.pi_red);  % Class probs 
    analysis.pi_med = analysis.pi_med / sum(analysis.pi_med);
    % Item-response probs 
    analysis.theta_med = reshape(median(analysis.theta_red), [data_vars.p, analysis.k_red, data_vars.d_max]);  
    analysis.theta_med = analysis.theta_med ./ sum(analysis.theta_med, 3);
    analysis.xi_med = median(analysis.xi_red);  % Probit coefs      
      
    % Update class membership probit component using unique classes and updated posterior median estimates for xi
    indiv_lik_probit_class = zeros(data_vars.n, analysis.k_red); % Initialize indiv probit lik matrix for each indiv and class
    analysis.z_i = MCMC_out.z_i(end, :)';                        % Obtain last values of probit latent variable
    for k = 1:analysis.k_red
        % Factor variable coding temporarily assign all indivs to class k
        Q_temp = zeros(data_vars.n, S * analysis.k_red);  % Initialize temp design matrix for class membership
        for s = 1:S
            Q_temp(:, ((s-1) * analysis.k_red) + k) = q_dem(:, s);         % Temporarily assign all indivs to class k
        end
%                % Reference cell coding temporarily assign all indivs to class k
%                class_temp = ones(data_vars.n, 1) * k;
%                q_design_temp = horzcat(q_dem, class_temp);
%                Q_temp = x2fx(q_design_temp, 'interaction', [1 2], [length(unique(q_dem)), analysis.k_red]);
%                q_class = zeros(data_vars.n, analysis.k_red);        % Initialize design matrix for class membership
%                q_class(:, k) = ones(data_vars.n, 1);                % Temporarily assign all indivs to class k
%                Q_temp = [q_dem q_class];                            % Form temporary covariate design matrix
       lin_pred_temp = Q_temp * transpose(analysis.xi_med); % Linear predictor, Q*xi
       % Update indiv probit likelihood assuming class k
       indiv_lik_probit_class(:, k) = normpdf(analysis.z_i,lin_pred_temp,1).*((data_vars.y==1).*(analysis.z_i>0)+(data_vars.y==0).*(analysis.z_i<=0));
    end  
    
    % Update class membership, c_i, using unique classes and updated posterior median estimates for pi and theta
 	prob_ci_numer = zeros(data_vars.n, analysis.k_red);  % Initialize indiv OFMM lik matrix for each indiv and class
    for k = 1:analysis.k_red
        % Matrix of unique item-level probs for class k
        theta_k = reshape(analysis.theta_med(:, k, :), data_vars.p, data_vars.d_max); 
        % Matrix of categ likelihood \prod_{r=1}^d \theta_{jr|k}^{I(x_ij=r)} for each indiv and item 
        categ_lik = reshape(theta_k(data_vars.lin_idx), [data_vars.n, data_vars.p]); 
        % Update numerator for posterior P(c_i=k|-) (indiv OFMM likelihood)
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
        % Factor variable coding design matrix with updated class memberships using median estimates 
        Q_med = zeros(data_vars.n, S * analysis.k_red);
        x_ci = dummyvar(analysis.c_i);
        for s = 1:S
            for k = 1:analysis.k_red
                Q_med(:, ((s-1) * analysis.k_red) + k) = q_dem(:, s) .* x_ci(:, k);
            end
        end
%             % Reference cell coding design matrix with demographic covariates and initial latent classes
%             q_design = horzcat(q_dem, analysis.c_i);
%             Q_med = x2fx(q_design, 'interaction', [1 2], [length(unique(q_dem)), analysis.k_red]);
%             Q_med = [q_dem dummyvar(analysis.c_i)];  % Design matrix with updated class memberships using median estimates 
    analysis.Phi_med = normcdf(Q_med * transpose(analysis.xi_med)); % Linear predictor, Q*xi. Probit model mean using median estimates
    analysis.Phi_med(analysis.Phi_med == 0) = 1e-10;  % Adjust extremes
    analysis.Phi_med(analysis.Phi_med == 1) = 1 - 1e-10;
    
%     % Update latent probit variable, z_i, using unique classes and updated posterior median estimates for xi
%     analysis.z_i = zeros(data_vars.n, 1);       % Initialize latent probit variable, z_i
%     % For cases, z_i ~ TruncNormal(mean=Q*xi, var=1, low=0, high=Inf)
%     analysis.z_i(data_vars.y == 1) = truncnormrnd(1, analysis.Phi_med(data_vars.y == 1), 1, 0, Inf); 
%     % For controls, z_i ~ TruncNormal(mean=Q*xi, var=1, low=-Inf, high=0)
%     analysis.z_i(data_vars.y == 0) = truncnormrnd(1, analysis.Phi_med(data_vars.y == 0), 1, -Inf, 0); 
%     if sum(analysis.z_i == Inf) > 0                  % Control extremes
%         analysis.z_i(analysis.z_i == Inf) = norminv(1 - 1e-10);
%     end 
%     if sum(analysis.z_i == -Inf) > 0
%         analysis.z_i(analysisz_i == -Inf) = norminv(1e-10);
%     end      
    
    % Update individual log-likelihood using unique classes and updated posterior median estimates 
    % Categorical component
    item_idx = repmat(1:data_vars.p, data_vars.n, 1);     % nxp matrix of item id's. Replicates 1:p row n times
    item_idx = item_idx(:);                               % Concat by col into (np)x1 vector. 1(x n),2(x n),...,p(x n)
    class_idx = repmat(analysis.c_i, 1, data_vars.p);     % Replicate col vector c_i p times to form nxp matrix of class assigns
    class_idx = class_idx(:);                             % Concat by col into (np)x1 vector 
    x_idx = data_vars.food(:);                            % Concat responses by col in (np)x1 vector 
    % lin_idx: (npx1) vector of linear indices indicating value of pxKxd theta matrix, for each item and indiv
    lin_idx = sub2ind([data_vars.p, analysis.k_red, data_vars.d_max], item_idx, class_idx, x_idx);
    % nxp matrix of theta values for each indiv and item \prod_{r=1}^d \theta_{jr|k}^{I(x_ij=r, c_i=k)}
    theta_indiv = reshape(analysis.theta_med(lin_idx), [data_vars.n, data_vars.p]);
    % Probit component 
    probit_lik = normpdf(analysis.z_i, analysis.Phi_med, 1) .* ((data_vars.y == 1).* (analysis.z_i > 0) + (data_vars.y == 0) .* (analysis.z_i <= 0));
    % Indiv complete log-likelihood
    if size(analysis.pi_med, 2) == 1  % Adjustment in case pi is 1x1
        analysis.loglik_med = log(analysis.pi_med(analysis.c_i) .* prod(theta_indiv, 2) .* probit_lik);
    else
        analysis.loglik_med = log(analysis.pi_med(analysis.c_i) * prod(theta_indiv, 2) * probit_lik);
    end

    % DIC is a metric for model goodness of fit. It consists of two terms:
        % 1) median(loglik): calculate log-lik at each MCMC iter then get the median
        % 2) loglik_med: get posterior estimates of each param, then use these to calculate the log-lik
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