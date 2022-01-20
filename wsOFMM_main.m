%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted Supervised OFMM Main %
% Programmer: SW                %
% Data: Simulations             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% Load simulated data
sim_n = bb;    % Simulation number. 'bb' created through 'job_array.sh' script
scenario = 1;  % Weighting scenario
samp_data = importdata(strcat('simdata_wsRPC_scen', num2str(scenario), '_iter',num2str(sim_n),'.mat'));

%% Get data variables
k_max = 50;    % Upper limit for number of classes
data_vars = wtd_get_data_vars(samp_data);

%% Initialize priors and variables for OFMM model
sp_k = k_max;                   % Denom constant to restrict alpha size and num classes for OFMM model
alpha = ones(1, k_max) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
eta = ones(1, data_vars.d_max); % Hyperparam for item response probs. 
OFMM_params = init_OFMM_params(data_vars, k_max, alpha, eta);

%% Initialize priors and variables for probit model
q_dem = [];                                 % Matrix of demographic covariates in cell-means format. Default is empty.
S = size(q_dem, 2);                         % Number of demographic covariates in the probit model
p_cov = k_max + S;                          % Number of covariates in probit model
mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVGamma(shape=5/2, scale=5/2)
Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
probit_params = init_probit_params(data_vars, k_max, q_dem, mu_0, Sig_0);

%% Run adaptive sampler to obtain number of classes
n_runs = 25000;  % Number of MCMC iterations
burn = 15000;    % Burn-in period
thin = 5;        % Thinning factor
[MCMC_out, ~, ~] = wtd_run_MCMC(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_max, q_dem, p_cov, alpha, eta, mu_0, Sig_0, S);
k_fixed = ceil(median(sum(MCMC_out.pi > 0.05, 2))); % Obtain fixed number of classes to use in the fixed sampler

%% Run fixed sampler to obtain posteriors and save output
% Initialize OFMM model using fixed number of classes
sp_k = k_fixed;                   % Denom constant to restrict alpha size and num classes for OFMM model
alpha = ones(1, k_fixed) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
OFMM_params = init_OFMM_params(data_vars, k_fixed, alpha, eta);

% Initialize probit model using fixed number of classes
p_cov = k_fixed + S;                        % Number of covariates in probit model
mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVGamma(shape=5/2, scale=5/2)
Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
probit_params = init_probit_params(data_vars, k_fixed, q_dem, mu_0, Sig_0);

% Run MCMC algorithm using fixed number of classes
[MCMC_out, OFMM_params, probit_params] = wtd_run_MCMC(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_fixed, q_dem, p_cov, alpha, eta, mu_0, Sig_0, S);
% Save MCMC output
save(strcat('wsOFMM_MCMC_scen', num2str(scenario), '_iter', num2str(sim_n)), 'MCMC_out');

%% Post-processing to recalibrate labels and remove extraneous empty classes
post_MCMC_out = post_process(MCMC_out, data_vars, S);

%% Obtain posterior estimates, reduce number of classes, analyze results, and save output
analysis = analyze_results(post_MCMC_out, data_vars, q_dem, S, p_cov);
% Save parameter estimates and analysis results
save(strcat('wsOFMM_Results_scen', num2str(scenario), '_iter', num2str(sim_n)), 'post_MCMC_out', 'analysis');



%% Miscellaneous additional code
%     % Create and save figures
%     figure; %check mixing of pi parameter
%     plot(MCMC_out.pi)
%     saveas(gcf,'wsOFMM_pis.png')
%     
%     figure; %save dendrogram of hierarchical clustering
%     dendrogram(post_MCMC_out.tree);
%     saveas(gcf,'wsOFMM_dendrogram.png')
%
%     % Compute the MSE comparing the estimated probit model to the true model 
%     y_mse = immse(analysis.Phi_med, samp_data.true_Phi)  % Phi_true must be pulled from simulated dataset
% 
%     % MCMC iterations for testing code
%     n_runs = 50;  % Number of MCMC iterations
%     burn = 30;    % Burn-in period
%     thin = 2;        % Thinning factor
%     sim_n = 1;    % Simulation iteration
    