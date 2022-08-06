%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted Supervised OFMM Main %
% Programmer: SW                %
% Data: Simulations             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% wsOFMM_main_latent takes in scenario, sim_n, and samp_n as command line 
% array index arguments. It reads in the simulated dataset corresponding to 
% the given array and scenario, runs the wsOFMM model with the latent 
% variables, and saves the MCMC output and posterior results as files.
% Inputs:
%   scenario: Command line argument indicating weighting scenario
%   sim_n: Command line argument indicating simulation array index
%   samp: Command line argument indicating sample index
% Outputs: saves MCMC output in a file named 'wsOFMM_MCMC...' and saves 
% posterior results in a file named 'wsOFMM_results...'

%% PARAMETER SETUP
scenario = 7; %% SRS with n=4000
sim_n = 1; % Iter = 1
samp_n = 5; % Samples 2, 5 tend to not converge
rng(sim_n, 'twister');  % set seed

%% Load simulated data and check if file already exists
% Input directory
in_dir = strcat(pwd, '/');
% Output directory 
out_dir = in_dir;   
if samp_n > 0   % If working with sample 
    samp_data = importdata(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n), '.mat'));
else            % If working with population
    samp_data = importdata(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '.mat'));
end


%% FIXED SAMPLER
% Get data variables
data_vars = wtd_get_data_vars_latent(samp_data);
k_fixed = 2;
eta = ones(1, data_vars.d_max); % Hyperparam for item response probs. 
q_dem = dummyvar(samp_data.true_Si);        % Matrix of demographic covariates in cell-means format. Default contains subpop. 
S = size(q_dem, 2);                         % Number of demographic covariates in the probit model

% Initialize OFMM model using fixed number of classes 
sp_k = k_fixed;                   % Denom constant to restrict alpha size and num classes for OFMM model
alpha = ones(1, k_fixed) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
OFMM_params = wtd_init_OFMM_params_latent(data_vars, k_fixed, alpha, eta);

% Initialize probit model using fixed number of classes
p_cov = k_fixed + S;                        % Number of covariates in probit model
mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVGamma(shape=5/2, scale=5/2)
Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
probit_params = init_probit_params_latent(data_vars, k_fixed, q_dem, mu_0, Sig_0, OFMM_params);

% Run MCMC algorithm using fixed number of classes
n_runs = 2500;  % Number of MCMC iterations
burn = 1500;    % Burn-in period
thin = 5;        % Thinning factor
[MCMC_out, ~, ~] = wtd_run_MCMC_latent(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_fixed, q_dem, p_cov, alpha, eta, mu_0, Sig_0, S);

% Post-processing to recalibrate labels and remove extraneous empty classes
post_MCMC_out = post_process_latent(MCMC_out, data_vars, S);

%Obtain posterior estimates, reduce number of classes, analyze results, and save output
analysis = analyze_results_latent(MCMC_out, post_MCMC_out, data_vars, q_dem, S, p_cov);

%% CHECK CONVERGENCE
figure; %check mixing of pi parameter
plot(MCMC_out.pi)

figure; %check ordered pi
plot(post_MCMC_out.pi) 

figure; % check ordered theta
plot(post_MCMC_out.theta(:,1,1,1))
hold on
plot(post_MCMC_out.theta(:,1,1,2))
hold off

figure; % check ordered theta
plot(post_MCMC_out.theta(:,1,2,1))
hold on
plot(post_MCMC_out.theta(:,1,2,2))
hold off

figure; % check ordered xi
plot(post_MCMC_out.xi)

disp(samp_data.true_xi);
disp(analysis.xi_med);
tabulate(samp_data.true_Phi)
tabulate(analysis.Phi_med)
disp(mean(abs(samp_data.true_xi - analysis.xi_med)));
disp(mean(abs(samp_data.true_Phi - analysis.Phi_med)));
disp(sum(samp_data.true_xi - analysis.xi_med)^2);
disp(immse(samp_data.true_Phi, analysis.Phi_med));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADAPTIVE AND FIXED SAMPLERS
%% PARAMETER SETUP
scenario = 6; %% Wtd with n=4000
sim_n = 1; % Iter = 1
samp_n = 2; % Samples 2, 5 tend to not converge
rng(sim_n, 'twister');  % set seed

%% Load simulated data and check if file already exists
% Input directory
in_dir = strcat(pwd, '/');
% Output directory 
out_dir = in_dir;   
if samp_n > 0   % If working with sample 
    samp_data = importdata(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n), '.mat'));
else            % If working with population
    samp_data = importdata(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '.mat'));
end

% Get data variables
data_vars = wtd_get_data_vars_latent(samp_data);
k_max = 50;    % Upper limit for number of classes

% Initialize priors and variables for OFMM model
sp_k = k_max;                   % Denom constant to restrict alpha size and num classes for OFMM model
alpha = ones(1, k_max) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
eta = ones(1, data_vars.d_max); % Hyperparam for item response probs. 
OFMM_params = wtd_init_OFMM_params_latent(data_vars, k_max, alpha, eta);

% Initialize priors and variables for probit model
q_dem = dummyvar(samp_data.true_Si);        % Matrix of demographic covariates in cell-means format. Default contains subpop.                               
S = size(q_dem, 2);                         % Number of demographic covariates in the probit model
p_cov = k_max + S;                          % Number of covariates in probit model
mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVGamma(shape=5/2, scale=5/2)
Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
probit_params = init_probit_params_latent(data_vars, k_max, q_dem, mu_0, Sig_0, OFMM_params);

% Run adaptive sampler to obtain number of classes
n_runs = 2500;  % Number of MCMC iterations
burn = 1500;    % Burn-in period
thin = 5;        % Thinning factor
[MCMC_out, ~, ~] = wtd_run_MCMC_latent(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_max, q_dem, p_cov, alpha, eta, mu_0, Sig_0, S);

% Post-processing for adaptive sampler
post_MCMC_out = post_process_latent(MCMC_out, data_vars, S);
% Array of posterior median item-response probs
theta_med_temp = reshape(median(post_MCMC_out.theta), [data_vars.p, post_MCMC_out.k_med, data_vars.d_max]);  
% Matrix of most likely consump level based on item-response probs, for each item and class across MCMC runs
[~, ind0] = max(theta_med_temp, [], 3); 
t_ind0 = transpose(ind0);
% Identify unique classes   
[~, ia, ~] = unique(t_ind0, 'rows');  % Vector of class indices corresponding to unique classes
k_fixed = length(ia);          % Number of unique classes
clear OFMM_params probit_params MCMC_out post_MCMC_out;  % Reduce memory burden

% Run fixed sampler to obtain posteriors and save output
% Initialize OFMM model using fixed number of classes  
sp_k = k_fixed;                   % Denom constant to restrict alpha size and num classes for OFMM model
alpha = ones(1, k_fixed) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
OFMM_params = wtd_init_OFMM_params_latent(data_vars, k_fixed, alpha, eta);

% Initialize probit model using fixed number of classes
p_cov = k_fixed + S;                        % Number of covariates in probit model
mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVGamma(shape=5/2, scale=5/2)
Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
probit_params = init_probit_params_latent(data_vars, k_fixed, q_dem, mu_0, Sig_0, OFMM_params);

% Run MCMC algorithm using fixed number of classes
[MCMC_out, ~, ~] = wtd_run_MCMC_latent(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_fixed, q_dem, p_cov, alpha, eta, mu_0, Sig_0, S);

% Post-processing to recalibrate labels and remove extraneous empty classes
post_MCMC_out = post_process_latent(MCMC_out, data_vars, S);

%Obtain posterior estimates, reduce number of classes, analyze results, and save output
analysis = analyze_results_latent(MCMC_out, post_MCMC_out, data_vars, q_dem, S, p_cov);

%% CHECK CONVERGENCE
figure; %check mixing of pi parameter
plot(MCMC_out.pi)

figure; %check ordered pi
plot(post_MCMC_out.pi) 

figure; % check ordered theta
plot(post_MCMC_out.theta(:,1,1,1))
hold on
plot(post_MCMC_out.theta(:,1,1,2))
hold off

figure; % check ordered xi
plot(post_MCMC_out.xi)

figure; % check ordered theta
plot(post_MCMC_out.theta(:,1,2,1))
hold on
plot(post_MCMC_out.theta(:,1,2,2))
hold off
