%%%%%%%%%%%%%%%%%%%%%%
% Weighted OFMM Main %
% Programmer: SW     %
% Data: Simulations  %
%%%%%%%%%%%%%%%%%%%%%%

% wOFMM_main_latent takes in scenario, sim_n, and samp_n as command line 
% array index arguments. It reads in the simulated dataset corresponding to 
% the given array and scenario, runs the wsOFMM model with the latent 
% variables, and saves the MCMC output and posterior results as files.
% Inputs:
%   scenario: Command line argument indicating weighting scenario
%   sim_n: Command line argument indicating simulation array index
%   samp: Command line argument indicating sample index
% Outputs: saves posterior results in a file named 'wsOFMM_results...'
function wOFMM_main_latent(scenario, sim_n, samp_n)
    % set seed
    rng(sim_n, 'twister');

    %% Load simulated data and check if file already exists
    % Input directory
    in_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Data/";   
    % Output directory 
    out_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Results/";  
%             in_dir = strcat(pwd, '/');
%             out_dir = in_dir; 
    if samp_n > 0   % If working with sample 
        samp_data = importdata(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n), '.mat'));
        already_done = isfile(strcat(out_dir, 'wOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n), '.mat'));
%             already_done = false;
    else            % If working with population
        samp_data = importdata(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '.mat'));
        already_done = isfile(strcat(out_dir, 'wOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n), '.mat'));
    end
    
    % Run the model if the results file does not already exist
    if already_done
        % If the results file already exists, print out statement
        disp(strcat('Scenario ', num2str(scenario), ' iter ', num2str(sim_n), ' samp ', num2str(samp_n), ' already exists.'));
    else
        
        tic
        %% Get data variables
        data_vars = wtd_get_data_vars_latent(samp_data);
        k_max = 30;    % Upper limit for number of classes
        
        %% Initialize priors and variables for OFMM model
        sp_k = k_max;                   % Denom constant to restrict alpha size and num classes for OFMM model
        alpha = ones(1, k_max) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
        eta = ones(1, data_vars.d_max); % Hyperparam for item response probs. 
        OFMM_params = wtd_init_OFMM_params_latent(data_vars, k_max, alpha, eta);

        %% Run adaptive sampler to obtain number of classes
        n_runs = 25000;  % Number of MCMC iterations
        burn = 15000;    % Burn-in period
%                 n_runs = 250;  % Number of MCMC iterations
%                 burn = 150;    % Burn-in period        
        thin = 5;        % Thinning factor
        [MCMC_out, ~] = wtd_run_MCMC_unsup_latent(data_vars, OFMM_params, n_runs, burn, thin, k_max, alpha, eta);
        
        
        %% Post-processing for adaptive sampler
        % Post-processing for adaptive sampler
        post_MCMC_out = post_process_unsup_latent(MCMC_out, data_vars);
        % Array of posterior median item-response probs
        theta_med_temp = reshape(median(post_MCMC_out.theta), [data_vars.p, post_MCMC_out.k_med, data_vars.d_max]);  
        % Matrix of most likely consump level based on item-response probs, for each item and class across MCMC runs
        [~, ind0] = max(theta_med_temp, [], 3); 
        t_ind0 = transpose(ind0);
        % Identify unique classes   
        [~, ia, ~] = unique(t_ind0, 'rows');  % Vector of class indices corresponding to unique classes
        k_fixed = length(ia);          % Number of unique classes
        clear OFMM_params MCMC_out post_MCMC_out;  % Reduce memory burden


        %% Run fixed sampler to obtain posteriors and save output
        % Initialize OFMM model using fixed number of classes  
        sp_k = k_fixed;                   % Denom constant to restrict alpha size and num classes for OFMM model
        alpha = ones(1, k_fixed) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
        OFMM_params = wtd_init_OFMM_params_latent(data_vars, k_fixed, alpha, eta);

        % Run MCMC algorithm using fixed number of classes
        [MCMC_out, ~] = wtd_run_MCMC_unsup_latent(data_vars, OFMM_params, n_runs, burn, thin, k_fixed, alpha, eta); 

        %% Post-processing to recalibrate labels and remove extraneous empty classes
        post_MCMC_out = post_process_unsup_latent(MCMC_out, data_vars);

        %% Obtain posterior estimates, reduce number of classes, analyze results, FIT PROBIT MODEL
        analysis = analyze_results_unsup_latent(post_MCMC_out, data_vars);   
        
        Q_ref = table(categorical(samp_data.true_Si), categorical(analysis.c_i), data_vars.y);
        Q_ref.Properties.VariableNames = ["S_i", "c_i", "y"];
        fit = fitglm(Q_ref, 'y ~ S_i*c_i', 'Distribution', 'binomial', 'link', 'probit');
        coefs = fit.Coefficients.Estimate;
        ci = coefCI(fit);
        % Convert to factor variable coding from reference cell coding
        S = length(unique(samp_data.true_Si));
        analysis.xi_med = zeros(1, S * analysis.k_red);
        for s = 1:S
            for c = 1:analysis.k_red
                analysis.xi_med((s-1)*analysis.k_red + c) = coefs(1) + (c ~= 1)*coefs(S+(c-1)) + (s ~= 1)*coefs(s) + (c ~= 1)*(s ~= 1)*coefs(S + (analysis.k_red-1) + (c-1));
                analysis.xi_med_lb((s-1)*analysis.k_red + c) = ci(1,1) + (c ~= 1)*ci(S+(c-1),1) + (s ~= 1)*ci(s,1) + (c ~= 1)*(s ~= 1)*ci(S + (analysis.k_red-1) + (c-1),1);
                analysis.xi_med_ub((s-1)*analysis.k_red + c) = ci(1,2) + (c ~= 1)*ci(S+(c-1),2) + (s ~= 1)*ci(s,2) + (c ~= 1)*(s ~= 1)*ci(S + (analysis.k_red-1) + (c-1),2);
            end
        end
        
        
        runtime = toc;
        
%                 figure; %check ordered pi
%                 plot(post_MCMC_out.pi) 
%                 
%                 figure; % check ordered theta
%                 ylim([0,1])
%                 plot(post_MCMC_out.theta(:,1,1,1))
%                 hold on
%                 plot(post_MCMC_out.theta(:,1,1,2))
%                 hold on                 
%                 plot(post_MCMC_out.theta(:,1,1,3))
%                 hold on    
%                 plot(post_MCMC_out.theta(:,1,1,4))
%                 hold off                    

        %% Save parameter estimates and analysis results
        if samp_n > 0  % If working with sample 
            save(strcat(out_dir, 'wOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n)), 'post_MCMC_out', 'analysis', 'runtime');
        else
            save(strcat(out_dir, 'wOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n)), 'post_MCMC_out', 'analysis', 'runtime');
        end 
   end  
end


    
