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
%   samp_n: Command line argument indicating sample index
% Outputs: saves MCMC output in a file named 'sOFMM_MCMC...' and saves 
% posterior results in a file named 'sOFMM_results...'
function sOFMM_main_latent(scenario, sim_n, samp_n)
    % set seed
    rng(sim_n, 'twister');

    %% Load simulated data and check if file already exists
    % Input directory
    in_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Toy_Example/";
    % Output directory 
    out_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Toy_Example/";     
%             in_dir = strcat(pwd, '/');
%             out_dir = in_dir;  
    if samp_n > 0   % If working with sample 
        samp_data = importdata(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n), '.mat'));
        already_done = isfile(strcat(out_dir, 'sOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n), '.mat'));
%             already_done = false;
    else            % If working with population
        samp_data = importdata(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '.mat'));
        already_done = isfile(strcat(out_dir, 'sOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n), '.mat'));
    end
    
    % Run the model if the results file does not already exist
    if already_done
        % If the results file already exists, print out statement
        disp(strcat('Scenario ', num2str(scenario), ' iter ', num2str(sim_n), ' samp ', num2str(samp_n), ' already exists.'));
    else
        %% Get data variables
        data_vars = wtd_get_data_vars_latent(samp_data);  % Need weights for post-hoc
        k_max = 50;    % Upper limit for number of classes

        %% Initialize priors and variables for OFMM model
        sp_k = k_max;                   % Denom constant to restrict alpha size and num classes for OFMM model
        alpha = ones(1, k_max) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
        eta = ones(1, data_vars.d_max); % Hyperparam for item response probs. 
        OFMM_params = init_OFMM_params_latent(data_vars, k_max, alpha, eta);

        %% Initialize priors and variables for probit model
        % Factor variable version
        q_dem = dummyvar(samp_data.true_Si);        % Matrix of demographic covariates in cell-means format. Default contains subpop. 
        S = size(q_dem, 2);                         % Number of demographic covariates in the probit model
        p_cov = S * k_max;                          % Number of covariates; from interactions
%                 % Reference cell version
%                 q_dem = samp_data.true_Si;
%                 S = length(unique(q_dem));
%                 p_cov = S + k_max-1 + (S-1)*(k_max-1);  % Intercept + S_dummies + C_dummies + interactions
%                 % Cell means version
%                 p_cov = k_max + S;                          % Number of covariates in probit model
%        clear samp_data;                            % Reduce memory burden
        mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
        Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVInvGamma(shape=5/2, scale=5/2)
        Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
        probit_params = init_probit_params_latent(data_vars, k_max, q_dem, mu_0, Sig_0, OFMM_params);

        %% Run adaptive sampler to obtain number of classes
        n_runs = 25000;  % Number of MCMC iterations
        burn = 15000;    % Burn-in period
%                 n_runs = 2500;  % Number of MCMC iterations
%                 burn = 1500;    % Burn-in period
        thin = 5;        % Thinning factor
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
        p_cov = S * k_fixed;                          % Number of covariates; from interactions
%                 % Reference cell version
%                 p_cov = S + k_fixed-1 + (S-1)*(k_fixed-1);  % Intercept + S_dummies + C_dummies + interactions
%                 % Cell means vesion
%                 p_cov = k_fixed + S;                          % Number of covariates in probit model
        mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
        Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVGamma(shape=5/2, scale=5/2)
        Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
        probit_params = init_probit_params_latent(data_vars, k_fixed, q_dem, mu_0, Sig_0, OFMM_params);

        % Run MCMC algorithm using fixed number of classes
        [MCMC_out, ~, ~] = run_MCMC_latent(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_fixed, q_dem, p_cov, alpha, eta, mu_0, Sig_0, S);
%         % Save output
%         if samp_n > 0  % If working with sample 
%             save(strcat(out_dir, 'sOFMM_latent_MCMC_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n)), 'MCMC_out');
%         else
%             save(strcat(out_dir, 'sOFMM_latent_MCMC_scen', num2str(scenario), '_iter', num2str(sim_n)), 'MCMC_out');
%         end

        %% Post-processing to recalibrate labels and remove extraneous empty classes
        post_MCMC_out = post_process_latent(MCMC_out, data_vars, S);

        %% Obtain posterior estimates, reduce number of classes, analyze results, and save output
        analysis = analyze_results_latent(MCMC_out, post_MCMC_out, data_vars, q_dem, S, p_cov);
        
        % Post-hoc pi correction
        analysis.posthoc_pi = zeros(analysis.k_red, 1);
        for k = 1:analysis.k_red  % For each latent class
            % Get updated weighted num indivs assigned to class k
            analysis.posthoc_pi(k) = sum(data_vars.wt_kappa(analysis.c_i == k)) / sum(data_vars.wt_kappa);  
        end
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
                
        % Save parameter estimates and analysis results
        if samp_n > 0  % If working with sample 
            save(strcat(out_dir, 'sOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n)), 'post_MCMC_out', 'analysis');
        else
            save(strcat(out_dir, 'sOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n)), 'post_MCMC_out', 'analysis');
        end    
    end    
end


%% Miscellaneous additional code
%     % Create and save figures
%    figure; %check mixing of pi parameter
%    plot(MCMC_out.pi)
%    saveas(gcf,'wsOFMM_pis.png')
%     
%    figure; %save dendrogram of hierarchical clustering
%    dendrogram(post_MCMC_out.tree)
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
    
% %% CHECK CONVERGENCE
% figure; %check mixing of pi parameter
% plot(MCMC_out.pi)
% 
% figure; %check ordered pi
% plot(post_MCMC_out.pi) 
% 
% figure; % check ordered theta
% plot(post_MCMC_out.theta(:,1,1,1))
% hold on
% plot(post_MCMC_out.theta(:,1,1,2))
% hold off
% 
% figure; % check ordered xi
% plot(post_MCMC_out.xi)
% 
% figure; % check ordered theta
% plot(post_MCMC_out.theta(:,1,2,1))
% hold on
% plot(post_MCMC_out.theta(:,1,2,2))
% hold off
% 
% immse(samp_data.true_Phi, analysis.Phi_med)
% sum(abs(sort(samp_data.true_Phi) - sort(analysis.Phi_med)))
% 
% %% CONVERTING CELL MEANS TO REFERENCE CELL
% if samp_data.true_Ci(1) ~= analysis.c_i(1)
%     beta_1 = analysis.xi_med(2) + analysis.xi_med(3);
%     beta_2 = analysis.xi_med(1) - analysis.xi_med(2);
%     beta_3 = analysis.xi_med(4) - analysis.xi_med(3);
%     beta_4 = analysis.xi_med(1) + analysis.xi_med(4) - beta_1 - beta_2 - beta_3;    
% else
%     beta_1 = analysis.xi_med(2) + analysis.xi_med(4);
%     beta_2 = analysis.xi_med(1) - analysis.xi_med(2);
%     beta_3 = analysis.xi_med(3) - analysis.xi_med(4);
%     beta_4 = analysis.xi_med(1) + analysis.xi_med(3) - beta_1 - beta_2 - beta_3;
% end
% beta = [beta_1 beta_2 beta_3 beta_4]
% % Find individuals with the four patterns
% Q1 = find((samp_data.true_Si == 2) & (samp_data.true_Ci == 2));
% Q2 = find((samp_data.true_Si == 1) & (samp_data.true_Ci == 2));
% Q3 = find((samp_data.true_Si == 2) & (samp_data.true_Ci == 1));
% Q4 = find((samp_data.true_Si == 1) & (samp_data.true_Ci == 1));
% Q_indivs = [Q1(1); Q2(1); Q3(1); Q4(1)];
% pred_Phi = analysis.Phi_med(Q_indivs);
% %pred_Phi = normcdf(Q_test_ref * transpose(beta));
% % Factor variable
% Q_test_ref = [1 0 0 0;   % S=1,C=1
%               0 1 0 0;   % S=1,C=2
%               0 0 1 0;   % S=2,C=1
%               0 0 0 1];  % S=2,C=2
% %     Q_test_ref = [1 0 0 0;   % S=2,C=2
% %                   1 1 0 0;   % S=1,C=2
% %                   1 0 1 0;   % S=2,C=1
% %                   1 1 1 1];  % S=1,C=1
% actual_Phi = normcdf(Q_test_ref * transpose(samp_data.true_xi));
% sum(abs(sort(actual_Phi) - sort(pred_Phi)))
% pred_Phi
% normcdf(Q_test_ref * transpose(beta))
% 
% sim_data = samp_data;
% % Obtain sensitivity and specificity
% pw_classes_sens = zeros(analysis.k_red, sim_data.true_K); % Initialize matrix of pairwise class assigns for sensitivity
% pw_classes_spec = zeros(analysis.k_red, sim_data.true_K); % Initialize matrix of pairwise class assigns for specificity
% for j = 1:sim_data.true_K                                 % For each true class
%     for i = 1:analysis.k_red                              % For each predicted class
%         % Obtain prop subjects in same class who are correctly predicted to share a class
%         pw_classes_sens(i, j) = sum(analysis.c_i == i & sim_data.true_Ci == j) / sum(sim_data.true_Ci == j);
%         % Obtain prop subjects NOT in same class who are correctly predicted to NOT share a class
%         pw_classes_spec(i, j) = sum(analysis.c_i ~= i & sim_data.true_Ci ~= j) / sum(sim_data.true_Ci ~= j);
%     end
% end
% sens_class = max(pw_classes_sens);  % Max per column gives best pairing of true and predicted class membership    
% spec_class = max(pw_classes_spec);  % Max per column gives best pairing of true and predicted class membership 
% sens = mean(sens_class) % Mean sensitivity over all true classes
% spec = mean(spec_class) % Mean specificity over all true classes
% 
% % Find incorrectly classified individuals
% if samp_data.true_Ci(1) ~= analysis.c_i(1)
%     temp_ci = ones(size(analysis.c_i));
%     temp_ci(analysis.c_i == 1) = 2;
%     mismatch = find(temp_ci ~= samp_data.true_Ci);
% else
%     mismatch = find(analysis.c_i ~= samp_data.true_Ci);
% end    
% length(mismatch)
% 
% 
% 
% % Find differences using observed Phi values
% true_Phi_values = zeros(length(samp_data.true_xi), 1);
% % P(Y=1|S=1,C=1)
% subset = samp_data.Y_data(samp_data.true_Si == 1 & samp_data.true_Ci == 1);
% true_Phi_values(1) = sum(subset) / length(subset);
% % P(Y=1|S=1,C=2)
% subset = samp_data.Y_data(samp_data.true_Si == 1 & samp_data.true_Ci == 2);
% true_Phi_values(2) = sum(subset) / length(subset);
% % P(Y=1|S=2,C=1)
% subset = samp_data.Y_data(samp_data.true_Si == 2 & samp_data.true_Ci == 1);
% true_Phi_values(3) = sum(subset) / length(subset);
% % P(Y=1|S=2,C=2)
% subset = samp_data.Y_data(samp_data.true_Si == 2 & samp_data.true_Ci == 2);
% true_Phi_values(4) = sum(subset) / length(subset);
% true_Phi_values
% sum(abs(sort(true_Phi_values) - sort(pred_Phi)))
% 
% % sum abs dev and sum squared distance for xi's
% sum(abs(sort(samp_data.true_xi) - sort(analysis.xi_med)))
% sum((sort(samp_data.true_xi) - sort(analysis.xi_med)).^2) % take this over iterations to get MSE
