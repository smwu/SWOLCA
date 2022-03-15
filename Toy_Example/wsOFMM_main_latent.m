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
function wsOFMM_main_latent(scenario, sim_n, samp_n)
    % set seed
    rng(sim_n, 'twister');

    %% Load simulated data and check if file already exists
    % Input directory
    in_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Toy_Example/";   
    % Output directory 
    out_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Toy_Example/";  
    if samp_n > 0   % If working with sample 
        samp_data = importdata(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n), '.mat'));
        already_done = isfile(strcat(out_dir, 'wsOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n), '.mat'));
    else            % If working with population
        samp_data = importdata(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '.mat'));
        already_done = isfile(strcat(out_dir, 'wsOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n), '.mat'));
    end
    
    % Run the model if the results file does not already exist
    if already_done
        % If the results file already exists, print out statement
        disp(strcat('Scenario ', num2str(scenario), ' iter ', num2str(sim_n), ' samp ', num2str(samp_n), ' already exists.'));
    else
        %% Get data variables
        k_max = 50;    % Upper limit for number of classes
        data_vars = wtd_get_data_vars_latent(samp_data);

        %% Initialize priors and variables for OFMM model
        sp_k = k_max;                   % Denom constant to restrict alpha size and num classes for OFMM model
        alpha = ones(1, k_max) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
        eta = ones(1, data_vars.d_max); % Hyperparam for item response probs. 
        OFMM_params = wtd_init_OFMM_params_latent(data_vars, k_max, alpha, eta);

        %% Initialize priors and variables for probit model
        q_dem = dummyvar(samp_data.true_Si);        % Matrix of demographic covariates in cell-means format. Default contains subpop.                               
        clear samp_data;                            % Reduce memory burden
        S = size(q_dem, 2);                         % Number of demographic covariates in the probit model
        p_cov = k_max + S;                          % Number of covariates in probit model
        mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
        Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVGamma(shape=5/2, scale=5/2)
        Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
        probit_params = init_probit_params_latent(data_vars, k_max, q_dem, mu_0, Sig_0, OFMM_params);

        %% Run adaptive sampler to obtain number of classes
        n_runs = 25000;  % Number of MCMC iterations
        burn = 15000;    % Burn-in period
        thin = 5;        % Thinning factor
        [MCMC_out, ~, ~] = wtd_run_MCMC_latent(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_max, q_dem, p_cov, alpha, eta, mu_0, Sig_0);
        k_fixed = round(median(sum(MCMC_out.pi > 0.05, 2))); % Obtain fixed number of classes to use in the fixed sampler
        clear OFMM_params probit_params MCMC_out;           % Reduce memory burden

        %% Run fixed sampler to obtain posteriors and save output
        % Initialize OFMM model using fixed number of classes  
        sp_k = k_fixed;                   % Denom constant to restrict alpha size and num classes for OFMM model
        alpha = ones(1, k_fixed) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
        OFMM_params = wtd_init_OFMM_params_latent(data_vars, k_fixed, alpha, eta);

        % Initialize probit model using fixed number of classes
        p_cov = k_fixed + S;                        % Number of covariates in probit model
        mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
        Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVInvGamma(shape=5/2, scale=5/2)
        Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
        probit_params = init_probit_params_latent(data_vars, k_fixed, q_dem, mu_0, Sig_0, OFMM_params);

        % Run MCMC algorithm using fixed number of classes
        [MCMC_out, ~, ~] = wtd_run_MCMC_latent(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_fixed, q_dem, p_cov, alpha, eta, mu_0, Sig_0);
        % Save MCMC output
        if samp_n > 0  % If working with sample 
            save(strcat(out_dir, 'wsOFMM_latent_MCMC_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n)), 'MCMC_out');
        else
            save(strcat(out_dir, 'wsOFMM_latent_MCMC_scen', num2str(scenario), '_iter', num2str(sim_n)), 'MCMC_out');
        end    

        %% Post-processing to recalibrate labels and remove extraneous empty classes
        post_MCMC_out = post_process_latent(MCMC_out, data_vars, S);
        clear OFMM_params probit_params;  % Reduce memory burden

        %% Obtain posterior estimates, reduce number of classes, analyze results, and save output
        analysis = analyze_results_latent(MCMC_out, post_MCMC_out, data_vars, q_dem, S, p_cov);
        % Save parameter estimates and analysis results
        if samp_n > 0  % If working with sample 
            save(strcat(out_dir, 'wsOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n)), 'post_MCMC_out', 'analysis');
        else
            save(strcat(out_dir, 'wsOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n)), 'post_MCMC_out', 'analysis');
        end 
   end  
end


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
%     % Testing code
%     in_dir = "C:/Users/Lang/Documents/Harvard/Research/Briana/supRPC/wsOFMM/Data/";
%     out_dir = "C:/Users/Lang/Documents/Harvard/Research/Briana/supRPC/wsOFMM/Results/";
%     n_runs = 50;  % Number of MCMC iterations
%     burn = 30;    % Burn-in period
%     thin = 2;        % Thinning factor
%     sim_n = 1;    % Simulation iteration
    
