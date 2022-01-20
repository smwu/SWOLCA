% init_probit_params initializes priors and variables for the probit model given hyperparameters
% Inputs: 
%   data_vars: output from wtd_get_data_vars or get_data_vars function
%   k_max: max number of latent classes
%   q_dem: matrix of demographic covariates to include in regression, in cell-means format
%   mu_0: vector of mean hyperparams for probit coefficients, xi; (p_cov)x1. 
%         Default is normrnd(0,1,[p_cov,1]). MVN(0,1) with indep components
%   Sig_0: matrix of variance hyperparams for probit coefficients, xi; (p_cov)x(p_cov).
%          Default is diag(1./gamrnd(5/2,2/5,[pcov,1])). Gamma(shape=5/2,scale=5/2) with indep components
% Outputs: probit_params structural array with the following fields:
%   xi: vector prior for probit coefficients; 1x(p_cov). Dist: MVN(mu_0, sig_0)
%   indiv_lik_probit: matrix of indiv probit lik components used for posterior class membership; nx(k_max)
function probit_params = init_probit_params(data_vars, k_max, q_dem, mu_0, Sig_0)
    % Prior for model params, drawn from MVN(mu0, Sig0)
    probit_params.xi = mvnrnd(mu_0, Sig_0);  
    
    % Probit likelihood component of posterior class membership P(c_i|-)
    probit_params.indiv_lik_probit = zeros(data_vars.n, k_max);    % Initialize matrix for individual probit likelihood contribution
    q_class = zeros(data_vars.n, k_max);                           % Initialize design matrix for class membership
    for k = 1:k_max
       q_class(:, k) = ones(data_vars.n, 1);                       % Temporarily assign all indivs to class k
       Q_temp = [q_dem q_class];                                   % Form temporary covariate design matrix
       Phi_temp = normcdf(Q_temp * transpose(probit_params.xi));   % Vector of P(y_i=1|Q)
       Phi_temp(Phi_temp == 1) = 1 - 1e-10;                        % Adjust extremes
       Phi_temp(Phi_temp == 0) = 1e-10;
       loglik_probit = data_vars.y .* log(Phi_temp) + (1-data_vars.y) .* log(1-Phi_temp); % Log-lik helps avoid rounding issues
       probit_params.indiv_lik_probit(:, k) = exp(loglik_probit);  % Probit likelihood for assigning indivs to class k
    end
end