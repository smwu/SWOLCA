% init_probit_params_latent initializes priors and variables for the probit 
% model given hyperparameters using the latent variable formulation
% Inputs: 
%   data_vars: output from wtd_get_data_vars or get_data_vars function
%   k_max: max number of latent classes
%   q_dem: matrix of demographic covariates to include in regression, in cell-means format
%   mu_0: vector of mean hyperparams for probit coefficients, xi; (p_cov)x1. 
%         Default is normrnd(0,1,[p_cov,1]). MVN(0,1) with indep components
%   Sig_0: matrix of variance hyperparams for probit coefficients, xi; (p_cov)x(p_cov).
%          Default is diag(1./gamrnd(5/2,2/5,[pcov,1])). Gamma(shape=5/2,scale=5/2) with indep components
%   OFMM_params: output from init_OFMM_params_latent function
% Outputs: probit_params structural array with the following fields:
%   xi: vector prior for probit coefficients; 1x(p_cov). Dist: MVN(mu_0, sig_0)
%   z_i: vector of latent variables in probit model; nx1
%   probit_lik: vector of indiv probit likelihoods with class assignments; nx1
%   indiv_lik_probit_class: matrix of indiv probit lik with assumed classes, used for posterior class membership; nx(k_max)
function probit_params = init_probit_params_latent(data_vars, k_max, q_dem, mu_0, Sig_0, OFMM_params)
    % Prior for model params, drawn from MVN(mu0, Sig0)
    probit_params.xi = mvnrnd(mu_0, Sig_0);  
    
    % Probit model latent variable z_i 
    probit_params.z_i = zeros(data_vars.n, 1);   % Initialize latent probit variable, z_i
    % Factor variable coding design matrix with demographic covariates and initial latent classes
    S = size(q_dem, 2); 
    Q = zeros(data_vars.n, S * k_max);
    for s = 1:S
        for k = 1:k_max
            Q(:, ((s-1) * k_max) + k) = q_dem(:, s) .* OFMM_params.x_ci(:, k);
        end
    end
%             % Reference cell coding design matrix with demographic covariates and initial latent classes
%             q_design = horzcat(q_dem, OFMM_params.c_i);
%             Q = x2fx(q_design, 'interaction', [1 2], [length(unique(q_dem)), k_max]);
%             Q = [q_dem OFMM_params.x_ci];                % Design matrix with demographic covariates and initial latent classes
    lin_pred = Q * transpose(probit_params.xi);  % Linear predictor, Q*xi. Mean of truncated normal dist
    % For cases, z_i ~ TruncNormal(mean=Q*xi, var=1, low=0, high=Inf)
    probit_params.z_i(data_vars.y == 1) = truncnormrnd(1, lin_pred(data_vars.y == 1), 1, 0, Inf); 
    % For controls, z_i ~ TruncNormal(mean=Q*xi, var=1, low=-Inf, high=0)
    probit_params.z_i(data_vars.y == 0) = truncnormrnd(1, lin_pred(data_vars.y == 0), 1, -Inf, 0); 
    if sum(probit_params.z_i == Inf) > 0                  % Control extremes
        probit_params.z_i(probit_params.z_i == Inf) = norminv(1 - 1e-10);
    end 
    if sum(probit_params.z_i == -Inf) > 0
        probit_params.z_i(probit_params.z_i == -Inf) = norminv(1e-10);
    end 
    
    % Individual probit log-likelihood with class assignments and latent variable formulation    
    probit_params.probit_lik = normpdf(probit_params.z_i, lin_pred, 1) .* ((data_vars.y == 1).* (probit_params.z_i > 0) + (data_vars.y == 0) .* (probit_params.z_i <= 0));
    
    % Probit likelihood component of posterior class membership P(c_i|-) with assumed classes
    probit_params.indiv_lik_probit_class = zeros(data_vars.n, k_max); % Initialize matrix of indiv probit likelihood assuming each class       
    for k = 1:k_max
        % Factor variable coding temporarily assign all indivs to class k
        Q_temp = zeros(data_vars.n, S * k_max);  % Initialize temp design matrix for class membership
        for s = 1:S
            Q_temp(:, ((s-1) * k_max) + k) = q_dem(:, s);         % Temporarily assign all indivs to class k
        end
%                % Reference cell coding temporarily assign all indivs to class k
%                class_temp = ones(data_vars.n, 1) * k;
%                q_design_temp = horzcat(q_dem, class_temp);
%                Q_temp = x2fx(q_design_temp, 'interaction', [1 2], [length(unique(q_dem)), k_max]);
%                q_class = zeros(data_vars.n, k_max);                  % Initialize design matrix for class membership
%                q_class(:, k) = ones(data_vars.n, 1);                 % Temporarily assign all indivs to class k
%                Q_temp = [q_dem q_class];                             % Form temporary covariate design matrix
       lin_pred_temp = Q_temp * transpose(probit_params.xi); % Linear predictor, Q*xi
       % Update indiv probit likelihood assuming class k
       probit_params.indiv_lik_probit_class(:, k) = normpdf(probit_params.z_i,lin_pred_temp,1).*((data_vars.y==1).*(probit_params.z_i>0)+(data_vars.y==0).*(probit_params.z_i<=0));
    end 
end