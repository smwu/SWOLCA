%%%%%%%%%%%%%%%%%%%%%%%%
% Supervised OFMM Main %
% Programmer: SW       %
% Data: Simulations    %
%%%%%%%%%%%%%%%%%%%%%%%%

% sOFMM_main_latent takes in scenario, sim_n, and samp_n as command line 
% array index arguments. It reads in the simulated dataset corresponding to 
% the given array and scenario, runs the sOFMM model with the latent 
% variables, and saves the MCMC output and posterior results as files.
% Inputs:
%   scenario: Command line argument indicating weighting scenario
%   sim_n: Command line argument indicating simulation array index
%   samp: Command line argument indicating sample index
% Outputs: saves MCMC output in a file named 'sOFMM_MCMC...' and saves 
% posterior results in a file named 'sOFMM_results...'

%% PARAMETER SETUP
scenario = 15; %% SRS with n=4000
sim_n = 1; % Iter = 1
samp_n = 5; % Samples 2, 5 tend to not converge
rng(sim_n, 'twister');  % set seed

scenario = 5; %% SRS with n=400
sim_n = 1; % Iter = 1
samp_n = 1; % Sample number = 1
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
% Reference cell version
q_dem = samp_data.true_Si;
S = length(unique(q_dem));

% Initialize OFMM model using fixed number of classes 
sp_k = k_fixed;                   % Denom constant to restrict alpha size and num classes for OFMM model
alpha = ones(1, k_fixed) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
OFMM_params = init_OFMM_params_latent(data_vars, k_fixed, alpha, eta);

% Initialize probit model using fixed number of classes
p_cov = k_fixed + S;                        % Number of covariates in probit model
% Reference cell version
p_cov = S + k_fixed-1 + (S-1)*(k_fixed-1);  % Intercept + S_dummies + C_dummies + interactions
mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVGamma(shape=5/2, scale=5/2)
Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
probit_params = init_probit_params_latent(data_vars, k_fixed, q_dem, mu_0, Sig_0, OFMM_params);

% Run MCMC algorithm using fixed number of classes
n_runs = 2500;  % Number of MCMC iterations
burn = 1500;    % Burn-in period
thin = 5;        % Thinning factor
[MCMC_out, ~, ~] = run_MCMC_latent(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_fixed, q_dem, p_cov, alpha, eta, mu_0, Sig_0, S);

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


Q_test_ref = [1 0 0 0;   % S=1,C=1
              1 1 0 0;   % S=2,C=1
              1 0 1 0;   % S=1,C=2
              1 1 1 1];  % S=2,C=2
pred_Phi = normcdf(Q_test_ref * transpose(analysis.xi_med));
actual_Phi = normcdf(Q_test_ref * transpose(samp_data.true_xi));
sum(abs(actual_Phi - pred_Phi))


immse(samp_data.true_Phi, analysis.Phi_med)
sum(abs(sim_data.true_Phi - analysis.Phi_med))




%% Get data variables
k_max = 50;    % Upper limit for number of classes
data_vars = wtd_get_data_vars_latent(samp_data);

%% Initialize priors and variables for OFMM model
sp_k = k_max;                   % Denom constant to restrict alpha size and num classes for OFMM model
alpha = ones(1, k_max) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
eta = ones(1, data_vars.d_max); % Hyperparam for item response probs. 
OFMM_params = init_OFMM_params_latent(data_vars, k_max, alpha, eta);

%% Initialize priors and variables for probit model
q_dem = dummyvar(samp_data.true_Si);        % Matrix of demographic covariates in cell-means format. Default contains subpop. 
clear samp_data;                            % Reduce memory burden
S = size(q_dem, 2);                         % Number of demographic covariates in the probit model
p_cov = k_max + S;                          % Number of covariates in probit model
mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVInvGamma(shape=5/2, scale=5/2)
Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
probit_params = init_probit_params_latent(data_vars, k_max, q_dem, mu_0, Sig_0, OFMM_params);

%% Run adaptive sampler to obtain number of classes
n_runs = 100000;  % Number of MCMC iterations
burn = 70000;    % Burn-in period
thin = 5;        % Thinning factor
[MCMC_out, ~, ~] = run_MCMC_latent(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_max, q_dem, p_cov, alpha, eta, mu_0, Sig_0, S);
k_fixed = round(median(sum(MCMC_out.pi > 0.05, 2))); % Obtain fixed number of classes to use in the fixed sampler
clear OFMM_params probit_params MCMC_out;           % Reduce memory burden

%% Run fixed sampler to obtain posteriors and save output
% Initialize OFMM model using fixed number of classes        
sp_k = k_fixed;                   % Denom constant to restrict alpha size and num classes for OFMM model
alpha = ones(1, k_fixed) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
OFMM_params = init_OFMM_params_latent(data_vars, k_fixed, alpha, eta);

% Initialize probit model using fixed number of classes
p_cov = k_fixed + S;                        % Number of covariates in probit model
mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVGamma(shape=5/2, scale=5/2)
Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
probit_params = init_probit_params_latent(data_vars, k_fixed, q_dem, mu_0, Sig_0, OFMM_params);

% Run MCMC algorithm using fixed number of classes
%         n_runs = 50000;  % Number of MCMC iterations
%         burn = 40000;    % Burn-in period
%         thin = 5;        % Thinning factor
[MCMC_out, ~, ~] = run_MCMC_latent(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_fixed, q_dem, p_cov, alpha, eta, mu_0, Sig_0, S);
%         % Save output
%         if samp_n > 0  % If working with sample 
%             save(strcat(out_dir, 'sOFMM_latent_MCMC_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n)), 'MCMC_out');
%         else
%             save(strcat(out_dir, 'sOFMM_latent_MCMC_scen', num2str(scenario), '_iter', num2str(sim_n)), 'MCMC_out');
%         end

%% Post-processing to recalibrate labels and remove extraneous empty classes
post_MCMC_out = post_process_latent(MCMC_out, data_vars, S);
clear OFMM_params probit_params;  % Reduce memory burden

%% Obtain posterior estimates, reduce number of classes, analyze results, and save output
analysis = analyze_results_latent(MCMC_out, post_MCMC_out, data_vars, q_dem, S, p_cov);

% Post-hoc pi correction
analysis.posthoc_pi = zeros(analysis.k_red, 1);
for k = 1:analysis.k_red  % For each latent class
    % Get updated weighted num indivs assigned to class k
    analysis.posthoc_pi(k) = sum(data_vars.wt_kappa(analysis.c_i == k)) / sum(data_vars.wt_kappa);  
end

% Save parameter estimates and analysis results
if samp_n > 0  % If working with sample 
    save(strcat(out_dir, 'sOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n)), 'post_MCMC_out', 'analysis');
else
    save(strcat(out_dir, 'sOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n)), 'post_MCMC_out', 'analysis');
end        


%% Miscellaneous additional code
%     % Create and save figures
    figure; %check mixing of pi parameter
    plot(MCMC_out.pi)
%    saveas(gcf,'wsOFMM_pis.png')
%     
    figure; %check mixing of theta parameter
    plot(MCMC_out.theta)

    figure; %save dendrogram of hierarchical clustering
    dendrogram(post_MCMC_out.tree)
%    saveas(gcf,'wsOFMM_dendrogram.png')
%
%     % Compute the MSE comparing the estimated probit model to the true model 
%     y_mse = immse(analysis.Phi_med, samp_data.true_Phi)  % Phi_true must be pulled from simulated dataset
% 
%     % Testing code
%     in_dir = strcat(pwd, "/Data/");
%     out_dir = strcat(pwd, "/Results/");
%     n_runs = 50;  % Number of MCMC iterations
%     burn = 30;    % Burn-in period
%     thin = 2;     % Thinning factor
%     sim_n = 1;    % Simulation iteration
    
