% wtd_run_MCMC_unsup_latent runs the Gibbs Sampler MCMC algorithm to obtain 
% posteriors for the wOFMM model
% Inputs: 
%   data_vars: output from wtd_get_data_vars_latent function
%   OFMM_params: output from init_OFMM_params_latent function
%   n_runs: number of MCMC iterations
%   burn: burn-in period
%   thin: thinning factor
%   k_max: max number of latent classes
%   alpha: vector hyperparam for class membership probs, pi; 1x(k_max)
%   eta: vector hyperparam for item-response probs, theta; 1x(d_max)
% Outputs: returns and saves MCMC_out structural array with the following fields:
%   pi: matrix of class membership probabilities; (n_runs/thin)x(k_max)  
%   theta: 4-D array of item-response probabilities; (n_runs/thin)xpx(k_max)x(d_max) 
%   c_i: matrix of class assignments; (n_runs/thin)xn
function [MCMC_out, OFMM_params] = wtd_run_MCMC_unsup_latent(data_vars, OFMM_params, n_runs, burn, thin, k_max, alpha, eta)

    % Initialize output variables
    MCMC_out.pi = zeros(n_runs / thin, k_max);  
    MCMC_out.theta = zeros(n_runs / thin, data_vars.p, k_max, data_vars.d_max); 
    MCMC_out.c_i = zeros(n_runs / thin, data_vars.n);  

    % Run MCMC iterations
    tic                  % Start timer
    for iter = 1:n_runs  % For each MCMC run
        % Update parameters and variables

        % Update posterior class membership probs, pi, incorporating normalized weights
        for k = 1:k_max  % For each latent class
            % Get updated weighted num indivs assigned to class k
            OFMM_params.n_ci(k) = sum(data_vars.wt_kappa(OFMM_params.c_i == k));  
        end
        alpha_post = alpha + OFMM_params.n_ci;    % Posterior hyperparam is alpha(k) + n_ci(k)
        OFMM_params.pi = drchrnd(alpha_post, 1);  % Draw from posterior dist: Dir(alpha_post)

        % Update class membership, c_i
        prob_ci_numer = zeros(data_vars.n, k_max);  % Initialize matrix of P(c_i=k|-) numerator for each indiv and class
        for k = 1:k_max
            % Matrix of unique item-level probs assuming class k
            theta_k = reshape(OFMM_params.theta(:, k, :), data_vars.p, data_vars.d_max); 
            % Matrix of categ likelihood \prod_{r=1}^d \theta_{jr|k}^{I(x_ij=r)|c_i=k} for each indiv and item assuming class k 
            categ_lik = reshape(theta_k(data_vars.lin_idx), [data_vars.n, data_vars.p]); 
            % Update numerator for posterior P(c_i=k|-) 
            prob_ci_numer(:, k) = OFMM_params.pi(k) * prod(categ_lik, 2);  
        end
        % Matrix of updated posterior P(c_i=k|-) for each indiv and class
        prob_ci = bsxfun(@times, prob_ci_numer, 1 ./ (sum(prob_ci_numer, 2)));
        OFMM_params.x_ci = mnrnd(1, prob_ci);  % Matrix of updated c_i drawn from Mult(1, prob_ci). Each row is draw for an indiv
        [row, col] = find(OFMM_params.x_ci);   % Row and col indices of each nonzero element of x_ci; col is class assignment for each indiv
        sorted = sortrows([row col], 1);       % Sort indices in ascending row index order, to match subject index
        OFMM_params.c_i = sorted(:, 2);        % Vector of updated class assignment, c_i, for each individual

        % Update posterior item response probabilities, theta, incorporating normalized weights    
        n_theta = zeros(data_vars.p, data_vars.d_max); % Initialize matrix of num indivs with each item response level
        for k = 1:k_max  
            c_i_mat = repmat(OFMM_params.c_i, [1, data_vars.p]);  % Rep c_i's p times to form matrix of classes for each indiv and item   
            class_wt = (c_i_mat == k) .* data_vars.wt_kappa_mat;  % Matrix of normalized weights for indivs assigned to class k. All other rows set to 0
            for r = 1:data_vars.d_max                             % For each consumption level r
                % Vector of num indivs with level r and class k, for each item (i.e., sum I(x_ij=r, c_i=k) over indivs)
                n_theta(:, r) = sum((data_vars.food == r) .* class_wt)'; 
            end
            for j = 1:data_vars.p                                     % For each food item
                dj = data_vars.d(j);                                  % Num response levels for item j
                eta_post = eta(1:dj) + n_theta(j, 1:dj);              % Posterior hyperparam is eta[r] + n_theta(j,r)
                OFMM_params.theta(j, k, 1:dj) = drchrnd(eta_post, 1); % Draw from posterior dist: Dir(eta), for every item and class
            end
        end           

        % Store posterior values if the iteration is selected based on thinning  
        if mod(iter, thin) == 0
            MCMC_out.pi(iter / thin, :) = OFMM_params.pi;
            MCMC_out.c_i(iter / thin, :) = OFMM_params.c_i;
            MCMC_out.theta(iter / thin, :, :, :) = OFMM_params.theta; 
        end


        % Relabel classes every 10 iterations to encourage mixing
        if mod(iter, 10) == 0
            new_order = randperm(k_max);                               % New labels for latent classes
            new_ci = OFMM_params.c_i;
            for k = 1:k_max
                new_ci(OFMM_params.c_i == k) = new_order(k);           % Class assignments with new labels
            end
            OFMM_params.c_i = new_ci;                               % Relabel class assignments
            OFMM_params.theta = OFMM_params.theta(:, new_order, :); % Relabel item-response probabilities
            OFMM_params.pi = OFMM_params.pi(new_order);                % Relabel class probabilities
        end
    end
    MCMC_out.runtime = toc;  % Stop timer

    % Discard burn-in
    MCMC_out.pi = MCMC_out.pi((burn / thin) + 1:end, :);
    MCMC_out.theta = MCMC_out.theta((burn / thin) + 1:end, :, :, :);
    MCMC_out.c_i = MCMC_out.c_i((burn / thin) + 1:end, :);

end