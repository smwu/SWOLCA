% wtd_run_MCMC runs the Gibbs Sampler MCMC algorithm to obtain posteriors
% for the wsOFMM model
% Inputs: 
%   data_vars: output from wtd_get_data_vars function
%   OFMM_params: output from init_OFMM_params function
%   probit_params: output from init_probit_params function
%   n_runs: number of MCMC iterations
%   burn: burn-in period
%   thin: thinning factor
%   k_max: max number of latent classes
%   q_dem: matrix of demographic covariates to include in regression, in cell-means format
%   p_cov: number of covariates in probit model
%   alpha: vector hyperparam for class membership probs, pi; 1x(k_max)
%   eta: vector hyperparam for item-response probs, theta; 1x(d_max)
%   mu_0: vector mean hyperparam for probit coefficients, xi; (p_cov)x1  
%   Sig_0: matrix variance hyperparam for probit coefficients, xi; (p_cov)x(p_cov)
%   S: number of demographic covariates in the probit model
% Outputs: returns and saves MCMC_out structural array with the following fields:
%   pi: matrix of class membership probabilities; (n_runs/thin)x(k_max)  
%   theta: 4-D array of item-response probabilities; (n_runs/thin)xpx(k_max)x(d_max) 
%   c_i: matrix of class assignments; (n_runs/thin)xn
%   xi: matrix of probit model coefficients; (n_runs/thin)x(p_cov) 
%   z_i: matrix of latent probit outcomes; (n_runs/thin)xn
%   loglik: vector of log likelihoods; (n_runs/thin)x1 
function [MCMC_out, OFMM_params, probit_params] = wtd_run_MCMC(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_max, q_dem, p_cov, alpha, eta, mu_0, Sig_0, S)
    % Initialize output variables
    MCMC_out.pi = zeros(n_runs / thin, k_max);  
    MCMC_out.theta = zeros(n_runs / thin, data_vars.p, k_max, data_vars.d_max); 
    MCMC_out.c_i = zeros(n_runs / thin, data_vars.n);  
    MCMC_out.xi = zeros(n_runs / thin, p_cov);  
    MCMC_out.z_i = zeros(n_runs / thin, data_vars.n);
    MCMC_out.loglik = zeros(n_runs / thin, 1); 
    
    % Run MCMC iterations
    tic                  % Start timer
    for iter = 1:n_runs  % For each MCMC run
        % Update parameters and variables
        [MCMC_out, OFMM_params, probit_params] = wtd_update_MCMC(MCMC_out, data_vars, OFMM_params, probit_params, thin, k_max, q_dem, alpha, eta, mu_0, Sig_0, iter);
    
        % Relabel classes every 10 iterations to encourage mixing
        if mod(iter, 10) == 0
            new_order = randperm(k_max);                               % New labels for latent classes
            new_ci = OFMM_params.c_i;
            for k = 1:k_max
                new_ci(OFMM_params.c_i == k) = new_order(k);           % Class assignments with new labels
            end
            OFMM_params.c_i = new_ci;                                  % Relabel class assignments
            OFMM_params.theta = OFMM_params.theta(:, new_order, :);    % Relabel item-response probabilities
            probit_params.indiv_lik_probit = probit_params.indiv_lik_probit(:, new_order);  % Relabel indiv probit likelihood
            % OFMM_params.pi = OFMM_params.pi(new_order);                % Relabel class probabilities
            % probit_params.xi((S+1):end) = probit_params.xi(S + new_order); % Relabel probit model coefficients
            % alpha = alpha(new_order);                                  % Relabel hyperparam for pi
            % mu_0((S+1):end) = mu_0(S + new_order);                     % Relabel mean hyperparam for xi
            % idx = find(Sig_0);                                         % Diagonal elements of Sig_0
            % Sig_0(idx((S+1):end)) = Sig_0(idx(S + new_order));         % Relabel var hyperparam for xi
        end
    end
    MCMC_out.runtime = toc;  % Stop timer
    
    % Discard burn-in
    MCMC_out.pi = MCMC_out.pi((burn / thin) + 1:end, :);
    MCMC_out.theta = MCMC_out.theta((burn / thin) + 1:end, :, :, :);
    MCMC_out.c_i = MCMC_out.c_i((burn / thin) + 1:end, :);
    MCMC_out.xi = MCMC_out.xi((burn / thin) + 1:end, :);  
    MCMC_out.z_i = MCMC_out.z_i((burn / thin) + 1:end, :);
    MCMC_out.loglik = MCMC_out.loglik((burn / thin) + 1:end); 
end