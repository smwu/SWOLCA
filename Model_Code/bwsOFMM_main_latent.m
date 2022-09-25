%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrap Weighted Supervised OFMM Main %
% Programmer: SW                          %
% Data: Simulations                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bwsOFMM_main_latent takes in scenario, sim_n, and samp_n as command line 
% array index arguments. It reads in the simulated dataset corresponding to 
% the given array and scenario, runs the wsOFMM model with the latent 
% variables, and saves the MCMC output and posterior results as files.
% Inputs:
%   scenario: Command line argument indicating weighting scenario
%   sim_n: Command line argument indicating simulation array index
%   samp_n: Command line argument indicating sample index
%   PRS_iter: Command line argument indicating PRS index
% Outputs: saves results in file named 'bwsOFMM_latent_results...'
function bwsOFMM_main_latent(scenario, sim_n, samp_n, PRS_iter)
    % set seed
    rng(sim_n + PRS_iter*1000, 'twister');

    %% Load simulated data and check if file already exists
    % Input directory
    in_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Data/";   
    % Output directory 
    out_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Results/M200/";  
%         in_dir = strcat(pwd, "/");
%         out_dir = strcat(pwd, "/");
    % Input directory
    if samp_n > 0   % If working with sample 
        samp_data = importdata(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n), '.mat'));
        already_done = isfile(strcat(out_dir, 'bwsOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n), '.mat'));
%             already_done = false;
    else            % If working with population
        samp_data = importdata(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '.mat'));
        already_done = isfile(strcat(out_dir, 'bwsOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n), '.mat'));
    end
    
    % Run the model if the results file does not already exist
    if already_done
        % If the results file already exists, print out statement
        disp(strcat('Scenario ', num2str(scenario), ' iter ', num2str(sim_n), ' samp ', num2str(samp_n), ' already exists.'));
    else
        
        tic
        %% Get data variables
        data_vars = get_data_vars_latent(samp_data);
        
        %% WFPBB
        N = sum(samp_data.sample_wt);  % Population size
        
        % Create synthetic population; Nxp
        synth_pop = wt_polya_post(samp_data.X_data, samp_data.sample_wt, data_vars.n, N, PRS_iter, samp_data.Y_data);
        % Form PRS by taking a simple random sample from the population
        prs_idx = datasample(1:N, data_vars.n, 'replace', false);
        
        % Set PRS as data to be used for the model
        samp_data.X_data = synth_pop.X_data(prs_idx, :);
        samp_data.Y_data = synth_pop.Y_data(prs_idx);
        data_vars = get_data_vars_latent(samp_data);
        
        %% Initialize priors and variables for OFMM model
        k_max = 30;    % Upper limit for number of classes
        sp_k = k_max;                   % Denom constant to restrict alpha size and num classes for OFMM model
        alpha = ones(1, k_max) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
        eta = ones(1, data_vars.d_max); % Hyperparam for item response probs. 
        OFMM_params = init_OFMM_params_latent(data_vars, k_max, alpha, eta);

        %% Initialize priors and variables for probit model
        % Factor variable version
        q_dem = dummyvar(samp_data.true_Si);        % Matrix of demographic covariates in cell-means format. Default contains subpop.                               
        S = size(q_dem, 2);                         % Number of demographic covariates in the probit model
        p_cov = k_max * S;                          % Number of covariates in probit model
        mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
        Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVGamma(shape=5/2, scale=5/2)
        Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
        probit_params = init_probit_params_latent(data_vars, k_max, q_dem, mu_0, Sig_0, OFMM_params);

        %% Run adaptive sampler to obtain number of classes
        n_runs = 25000;  % Number of MCMC iterations
        burn = 15000;    % Burn-in period
        thin = 5;        % Thinning factor
%             n_runs = 250;
%             burn = 150;
        [MCMC_out, ~, ~] = run_MCMC_latent(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_max, q_dem, p_cov, alpha, eta, mu_0, Sig_0, S);
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
                %k_fixed = round(median(sum(MCMC_out.pi > 0.05, 2))); % Obtain fixed number of classes to use in the fixed sampler
                %clear OFMM_params probit_params MCMC_out;           % Reduce memory burden
                
        %% Run fixed sampler to obtain posteriors and save output
        % Initialize OFMM model using fixed number of classes  
        sp_k = k_fixed;                   % Denom constant to restrict alpha size and num classes for OFMM model
        alpha = ones(1, k_fixed) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
        OFMM_params = init_OFMM_params_latent(data_vars, k_fixed, alpha, eta);

        % Initialize probit model using fixed number of classes
        % Factor variable version
        p_cov = S * k_fixed;                        % Number of covariates in probit model
        mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
        Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVInvGamma(shape=5/2, scale=5/2)
        Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
        probit_params = init_probit_params_latent(data_vars, k_fixed, q_dem, mu_0, Sig_0, OFMM_params);

        % Run MCMC algorithm using fixed number of classes
        [MCMC_out, ~, ~] = run_MCMC_latent(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_fixed, q_dem, p_cov, alpha, eta, mu_0, Sig_0, S);
%         % Save MCMC output
%         if samp_n > 0  % If working with sample 
%             save(strcat(out_dir, 'wsOFMM_latent_MCMC_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n)), 'MCMC_out');
%         else
%             save(strcat(out_dir, 'wsOFMM_latent_MCMC_scen', num2str(scenario), '_iter', num2str(sim_n)), 'MCMC_out');
%         end    

        %% Post-processing to recalibrate labels and remove extraneous empty classes
        post_MCMC_out = post_process_latent(MCMC_out, data_vars, S);

        %% Obtain posterior estimates, reduce number of classes, analyze results, and save output
        analysis = analyze_results_latent(MCMC_out, post_MCMC_out, data_vars, q_dem, S, p_cov);
        
        runtime = toc;
        
%                 figure; %check ordered pi
%                 plot(post_MCMC_out.pi) 
%                 
%                 figure; % check ordered theta
%                 plot(post_MCMC_out.theta(:,1,1,1))
%                 hold on
%                 plot(post_MCMC_out.theta(:,1,1,2))
%                 hold off
%                 
%                 figure; % check ordered xi
%                 plot(post_MCMC_out.xi)
% 
%                 sum(abs(sort(samp_data.true_xi) - sort(analysis.xi_med)))
%                 sum((sort(samp_data.true_xi) - sort(analysis.xi_med)).^2)

        % Save parameter estimates and analysis results
        if samp_n > 0  % If working with sample 
            save(strcat(out_dir, 'bwsOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n), '_prs', num2str(PRS_iter)), 'post_MCMC_out', 'analysis', 'runtime');
        else
            save(strcat(out_dir, 'bwsOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n), '_prs', num2str(PRS_iter)), 'post_MCMC_out', 'analysis', 'runtime');
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
    
