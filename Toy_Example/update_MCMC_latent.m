% update_MCMC_latent updates the posterior distributions of the parameters and 
% variables for the sOFMM model
% Inputs: 
%   MCMC_out: output initialized in run_MCMC function
%   data_vars: output from get_data_vars function
%   OFMM_params: output from init_OFMM_params function
%   probit_params: output from init_probit_params function
%   thin: thinning factor
%   k_max: max number of latent classes
%   q_dem: matrix of demographic covariates to include in regression, in cell-means format
%   alpha: vector hyperparam for class membership probs, pi
%   eta: vector hyperparam for item-response probs, theta
%   mu_0: vector mean hyperparam for probit coefficients, xi 
%   Sig_0: matrix variance hyperparam for probit coefficients, xi
%   iter: MCMC iteration indicator
% Outputs: returns updated OFMM_params and probit_params structural arrays,
% as well as updated MCMC_out structural array with the following fields:
%   pi: matrix of class membership probabilities 
%   theta: 4-D array of item-response probabilities
%   c_i: matrix of class assignments
%   xi: matrix of probit model coefficients
%   z_i: matrix of latent probit outcomes
%   loglik: vector of log likelihoods
function [MCMC_out, OFMM_params, probit_params] = update_MCMC_latent(MCMC_out, data_vars, OFMM_params, probit_params, thin, k_max, q_dem, alpha, eta, mu_0, Sig_0, iter) 
    % Update posterior class membership probs, pi
    for k = 1:k_max  % For each latent class
        % Get updated num indivs assigned to class k
        OFMM_params.n_ci(k) = sum(OFMM_params.c_i == k);  
    end
    alpha_post = alpha + OFMM_params.n_ci;    % Posterior hyperparam is alpha(k) + n_ci(k)
    OFMM_params.pi = drchrnd(alpha_post, 1);  % Draw from posterior dist: Dir(alpha_post) 
    
    
    
    % Update class membership, c_i
    prob_ci_numer_log = zeros(data_vars.n, k_max);  % Initialize log matrix to address underflow
 	%prob_ci_numer = zeros(data_vars.n, k_max);  % Initialize matrix of P(c_i=k|-) numerator for each indiv and class
    for k = 1:k_max
        % Matrix of unique item-level probs for class k
        theta_k = reshape(OFMM_params.theta(:, k, :), data_vars.p, data_vars.d_max); 
        % Matrix of categ likelihood \prod_{r=1}^d \theta_{jr|k}^{I(x_ij=r)} for each indiv and item 
        categ_lik = reshape(theta_k(data_vars.lin_idx), [data_vars.n, data_vars.p]); 
        % Update numerator for posterior P(c_i=k|-) 
        prob_ci_numer_log(:, k) = log(OFMM_params.pi(k)) + log(prod(categ_lik, 2)) + log(probit_params.indiv_lik_probit_class(:, k));
        %prob_ci_numer(:, k) = OFMM_params.pi(k) * prod(categ_lik, 2) .* probit_params.indiv_lik_probit_class(:, k);  
    end
    % Matrix of updated posterior P(c_i=k|-) for each indiv and class
    max_elems = max(prob_ci_numer_log, [], 2);
    subtract = prob_ci_numer_log - max_elems;
    log_denom = max_elems + log(sum(exp(subtract), 2));
    prob_ci = exp(prob_ci_numer_log - log_denom);
    %prob_ci = bsxfun(@times, prob_ci_numer, 1 ./ (sum(prob_ci_numer, 2)));
    % Adjust precision for mnrnd function
    
    
    
    
    
    
    
    % Update class membership, c_i
 	prob_ci_numer = zeros(data_vars.n, k_max);  % Initialize matrix of P(c_i=k|-) numerator for each indiv and class
    for k = 1:k_max
        % Matrix of unique item-level probs for class k
        theta_k = reshape(OFMM_params.theta(:, k, :), data_vars.p, data_vars.d_max); 
        % Matrix of categ likelihood \prod_{r=1}^d \theta_{jr|k}^{I(x_ij=r)} for each indiv and item 
        categ_lik = reshape(theta_k(data_vars.lin_idx), [data_vars.n, data_vars.p]); 
        % Update numerator for posterior P(c_i=k|-) 
        prob_ci_numer(:, k) = OFMM_params.pi(k) * prod(categ_lik, 2) .* probit_params.indiv_lik_probit_class(:, k);  
    end
    % Matrix of updated posterior P(c_i=k|-) for each indiv and class
    prob_ci = bsxfun(@times, prob_ci_numer, 1 ./ (sum(prob_ci_numer, 2)));
    % Adjust precision for mnrnd function
    for i = 1:data_vars.n
       prob_i = prob_ci(i, :);
       if any(prob_i < 1e-15)
           prob_i = prob_i + 1e-15;
           prob_ci(i, :) = prob_i ./ sum(prob_i);
       end    
    end       
    OFMM_params.x_ci = mnrnd(1, prob_ci);                   % Matrix of updated c_i drawn from Mult(1, prob_ci). Each row is draw for an indiv
    [row, col] = find(OFMM_params.x_ci == 1);               % Row and col indices of each nonzero element of x_ci; col is class assignment for each indiv
    sorted = sortrows([row col], 1);            % Sort indices in ascending row index order, to match subject index
    OFMM_params.c_i = sorted(:, 2);             % Vector of updated class assignment, c_i, for each individual
      
    % Update posterior item response probabilities, theta   
    n_theta = zeros(data_vars.p, data_vars.d_max); % Initialize matrix of num indivs with each item response level
    for k = 1:k_max  
        c_i_mat = repmat(OFMM_params.c_i, [1, data_vars.p]);  % Rep c_i's p times to form matrix of classes for each indiv and item
        class_wt = (c_i_mat == k);                            % Matrix indicating indivs assigned to class k. All other rows set to 0
        for r = 1:data_vars.d_max                             % For each consumption level r
            
            % DEBUGGING ERROR
            if any(size((data_vars.food == r)) ~= size(class_wt))
                disp("FOOD");
                disp(size((data_vars.food == r)));
                disp("CLASS_WT");
                disp(size(class_wt));
                disp("SORTED");
                disp(size(sorted));
                disp(sorted);
            end
            
            % Vector of num indivs with level r and class k, for each item (i.e., sum I(x_ij=r, c_i=k) over indivs)
            n_theta(:, r) = sum((data_vars.food == r) .* class_wt)'; 
        end
        for j = 1:data_vars.p                                     % For each food item
            dj = data_vars.d(j);                                  % Num response levels for item j
            eta_post = eta(1:dj) + n_theta(j, 1:dj);              % Posterior hyperparam is eta[r] + n_theta(j,r)
            OFMM_params.theta(j, k, 1:dj) = drchrnd(eta_post, 1); % Draw from posterior dist: Dir(eta), for every item and class
        end
    end        
    
    % Update latent probit variable, z_i
    Q = [q_dem OFMM_params.x_ci];                            % Design matrix with demographic covariates and updated latent classes
    lin_pred = Q * transpose(probit_params.xi);  % Linear predictor, Q*xi. Mean of truncated normal dist
    % For cases, z_i ~ TruncNormal(mean=Q*xi, var=1, low=0, high=Inf)
    probit_params.z_i(data_vars.y == 1) = truncnormrnd(1, lin_pred(data_vars.y == 1), 1, 0, Inf); 
    % For controls, z_i ~ TruncNormal(mean=Q*xi, var=1, low=-Inf, high=0)
    probit_params.z_i(data_vars.y == 0) = truncnormrnd(1, lin_pred(data_vars.y == 0), 1, -Inf, 0); 
    if sum(probit_params.z_i == Inf) > 0                       % Control extremes
        probit_params.z_i(probit_params.z_i == Inf) = norminv(1 - 1e-10);
    end 
    if sum(probit_params.z_i == -Inf) > 0
        probit_params.z_i(probit_params.z_i == -Inf) = norminv(1e-10);
    end    
    probit_params.probit_lik = normpdf(probit_params.z_i, lin_pred, 1) .* ((data_vars.y == 1).* (probit_params.z_i > 0) + (data_vars.y == 0) .* (probit_params.z_i <= 0));

    % Update individual log-likelihood with class assignments and latent variable formulation       
    item_idx = repmat(1:data_vars.p, data_vars.n, 1);     % nxp matrix of item id's. Replicates 1:p row n times
    item_idx = item_idx(:);                               % Concat by col into (np)x1 vector. 1(x n),2(x n),...,p(x n)
    class_idx = repmat(OFMM_params.c_i, 1, data_vars.p);  % Replicate col vector c_i p times to form nxp matrix of class assigns
    class_idx = class_idx(:);                             % Concat by col into (np)x1 vector 
    x_idx = data_vars.food(:);                            % Concat responses by col in (np)x1 vector 
    % lin_idx: (npx1) vector of linear indices indicating value of pxKxd theta matrix, for each item and indiv
    lin_idx = sub2ind([data_vars.p, k_max, data_vars.d_max], item_idx, class_idx, x_idx);
    % nxp matrix of theta values for each indiv and item \prod_{r=1}^d \theta_{jr|k}^{I(x_ij=r, c_i=k)}
    theta_indiv = reshape(OFMM_params.theta(lin_idx), [data_vars.n, data_vars.p]);
    % Indiv complete log-likelihood
    indiv_loglik = log(OFMM_params.pi(OFMM_params.c_i) * prod(theta_indiv, 2) * probit_params.probit_lik);
    
    % Update posterior probit model coefficients, xi   
    Sig_post = inv(Sig_0) + (transpose(Q) * Q);           % Precision of posterior Normal dist for xi
    mu_right = (Sig_0 \ mu_0) + (transpose(Q) * probit_params.z_i);       % Right part of mean of posterior Normal dist for xi
    % mu_right = (Sig_0\probit_params.xi0') + (transpose(Q)*z_i);  CHANGE THIS
    mu_post = Sig_post \ mu_right;                      % Mean of posterior Normal dist for xi
    probit_params.xi = mvnrnd(mu_post, inv(Sig_post));  % Draw from posterior dist: N(mu_post, Sig_post^{-1})    
    
    % Update probit likelihood component of posterior class membership P(c_i|-)
    for k = 1:k_max
       q_class = zeros(data_vars.n, k_max);                  % Initialize design matrix for class membership
       q_class(:, k) = ones(data_vars.n, 1);                 % Temporarily assign all indivs to class k
       Q_temp = [q_dem q_class];                             % Form temporary covariate design matrix
       lin_pred_temp = Q_temp * transpose(probit_params.xi); % Linear predictor, Q*xi       
       % Update indiv probit likelihood assuming class k
       probit_params.indiv_lik_probit_class(:, k) = normpdf(probit_params.z_i,lin_pred_temp,1).*((data_vars.y==1).*(probit_params.z_i>0)+(data_vars.y==0).*(probit_params.z_i<=0));
    end    
    
    % Store posterior values if the iteration is selected based on thinning  
    if mod(iter, thin) == 0
        MCMC_out.pi(iter / thin, :) = OFMM_params.pi;
        MCMC_out.c_i(iter / thin, :) = OFMM_params.c_i;
        MCMC_out.theta(iter / thin, :, :, :) = OFMM_params.theta; 
        MCMC_out.loglik(iter / thin) = sum(indiv_loglik); 
        MCMC_out.z_i(iter / thin, :) = probit_params.z_i;
        MCMC_out.xi(iter / thin, :) = probit_params.xi; 
    end
end
