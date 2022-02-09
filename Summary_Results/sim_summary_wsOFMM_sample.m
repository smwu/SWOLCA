%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Summary Results for sOFMM and wsOFMM
% when Dataset is a Sample
% Programmer: SW                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs:
%   scenario: Simulation scenario
%   iter: Vector of simulation iterations to aggregate over
%   samp_n: Vector of sample iterations to aggregate over. One of 'iter' and 'samp_n' must be scalar.
%   Model: Model to fit; either "wsOFMM" or "sOFMM"
function sim_summary_wsOFMM_sample(scenario, iter, samp_n, model)
    % Input directory for data
    data_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Data/";   
    % Input directory for model results 
    model_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Results/";      
    % Output directory 
    out_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Summary_Results/";  
    
    if length(samp_n) > length(iter)  % If aggregating over samples instead of iterations
        reps = samp_n;                % Replicates defined by sample iterations
    else                              % Aggregating over iterations
        reps = iter;                  % Replicates defined by simulation iterations
    end
    n_reps = length(reps);            % Number of replicates to aggregate across
    
    % Initialize outputs, stored for all simulated sets
    all.K = NaN(n_reps, 1);         % Initialize estimated number of classes
    all.K_match = NaN(n_reps, 1);   % Initialize boolean vector of correct estimation of number of classes
    all.sens = NaN(n_reps, 1);      % Initialize sensitivity
    all.spec = NaN(n_reps, 1);      % Initialize specificity
    all.dic6 = NaN(n_reps, 1);      % Initialize DIC version 6 to penalize overfitting 
    all.aebic = NaN(n_reps, 1);     % Initialize AEBIC model selection metric
    all.mse_Phi = NaN(n_reps, 1);   % Initialize MSE of class-specific response proportions
    all.mse_theta = NaN(n_reps, 1); % Initialize MSE of class-specific item response probs
    all.corr_S_C = NaN(n_reps, 1);  % Initialize correlation between S and C
    all.corr_S_Y = NaN(n_reps, 1);  % Initialize correlation between S and Y
    all.corr_C_Y = NaN(n_reps, 1);  % Initialize correlation between C and Y

    for i = 1:n_reps      % For each simulated set or sample
        sim_n = reps(i);  % Get iteration or sample index
        
        % Construct file names
        if length(samp_n) > length(iter)  % If aggregating over samples instead of iterations
            file_str = strcat(num2str(scenario), '_iter', num2str(iter), '_samp', num2str(samp_n(sim_n)), '.mat');
        else                              % Aggregating over iterations
            file_str = strcat(num2str(scenario), '_iter', num2str(iter(sim_n)), '_samp', num2str(samp_n), '.mat');
        end 
        
        % Check if the results file exists
        if isfile(strcat(model_dir, model, '_results_scen', file_str))

            % Load simulated dataset
            load(strcat(data_dir, 'simdata_scen', file_str), 'sim_data');
            % Load model results
            load(strcat(model_dir, model, '_results_scen', file_str), 'analysis') 

            all.K(sim_n) = analysis.k_red;                          % Estimated number of classes
            all.K_match(sim_n) = analysis.k_red == sim_data.true_K; % Estimated and true number of classes match (1/0)

            % Obtain sensitivity and specificity
            pw_classes_sens = zeros(analysis.k_red, sim_data.true_K); % Initialize matrix of pairwise class assigns for sensitivity
            pw_classes_spec = zeros(analysis.k_red, sim_data.true_K); % Initialize matrix of pairwise class assigns for specificity
            for j = 1:sim_data.true_K                                 % For each true class
                for i = 1:analysis.k_red                              % For each predicted class
                    % Obtain prop subjects in same class who are correctly predicted to share a class
                    pw_classes_sens(i, j) = sum(analysis.c_i == i & sim_data.true_Ci == j) / sum(sim_data.true_Ci == j);
                    % Obtain prop subjects NOT in same class who are correctly predicted to NOT share a class
                    pw_classes_spec(i, j) = sum(analysis.c_i ~= i & sim_data.true_Ci ~= j) / sum(sim_data.true_Ci ~= j);
                end
            end
            sens_class = max(pw_classes_sens);  % Max per column gives best pairing of true and predicted class membership    
            spec_class = max(pw_classes_spec);  % Max per column gives best pairing of true and predicted class membership 
            all.sens(sim_n) = mean(sens_class); % Mean sensitivity over all true classes
            all.spec(sim_n) = mean(spec_class); % Mean specificity over all true classes

            all.dic6(sim_n) = analysis.dic6;    % DIC version 6 to penalize overfitting
            all.aebic(sim_n) = analysis.aebic;  % AEBIC model selection metric

            % MSE of class-specific response proportions
            all.mse_Phi(sim_n) = immse(analysis.Phi_med, sim_data.true_Phi);  

            % MSE of class-specific item response probabilities
            pw_classes_mse = zeros(analysis.k_red, sim_data.true_K);         % Initialize matrix of pairwise class MSEs
            dim = size(sim_data.true_global_thetas); p = dim(1); d = dim(3); % Get food data dimensions
            for j = 1:sim_data.true_K                                        % For each true class
                % Obtain the true item response probs for true class j
                global_thetas_j = reshape(sim_data.true_global_thetas(:, j, :), [p, d]); 
                for i = 1:analysis.k_red                                     % For each predicted class
                    % Obtain the predicted item response probs for predicted class i
                    theta_med_i = reshape(analysis.theta_med(:, i, :), [p, d]);
                    % Obtain MSE of true and predicted item response probs for each pairwise class combo
                    pw_classes_mse(i, j) = immse(theta_med_i, global_thetas_j);
                end
            end
            mse_theta_class = min(pw_classes_mse);        % Min per column gives best pairing of true and predicted class membership    
            all.mse_theta(sim_n) = sum(mse_theta_class); % Total MSE is sum over all classes   

            % Get simulated correlations
            all.corr_S_C(sim_n) = corr(sim_data.true_Si, sim_data.true_Ci);
            all.corr_S_Y(sim_n) = corr(sim_data.true_Si, sim_data.Y_data);
            all.corr_C_Y(sim_n) = corr(sim_data.true_Ci, sim_data.Y_data);

        else
            % If the results file does not exist, print out statement
            disp(strcat('Scenario ', num2str(scenario), ' iter ', num2str(iter), ' samp ', num2str(sim_n), ' is missing.'));
        end      
    end
    
    % Median and 95% CI of estimated number of classes over all simulated datasets
    res.K_median = [median(all.K, 'omitnan') quantile(all.K, 0.025) quantile(all.K, 0.975)];          
    % Proportion of simulations with number of classes correctly estimated
    res.K_match = mean(all.K_match, 'omitnan');                     
    % Median and 95% CI of sensitivity over all simulated datasets
    res.sens_median = [median(all.sens, 'omitnan') quantile(all.sens, 0.025) quantile(all.sens, 0.975)];  
    % Median and 95% CI of specificity over all simulated datasets
    res.spec_median = [median(all.spec, 'omitnan') quantile(all.spec, 0.025) quantile(all.spec, 0.975)];
    % Median and 95% CI of DIC-6 over all simulated datasets
    res.dic6_median = [median(all.dic6, 'omitnan') quantile(all.dic6, 0.025) quantile(all.dic6, 0.975)];   
    % Median and 95% CI of AEBIC over all simulated datasets
    res.aebic_median = [median(all.aebic, 'omitnan') quantile(all.aebic, 0.025) quantile(all.aebic, 0.975)];  
    % Median and 95% CI of MSE of class-specific response proportions over all simulated datasets
    res.mse_Phi_median = [median(all.mse_Phi, 'omitnan') quantile(all.mse_Phi, 0.025) quantile(all.mse_Phi, 0.975)]; 
    % Median and 95% CI of MSE of item response probabilities over all simulated datasets 
    res.mse_theta_median = [median(all.mse_theta, 'omitnan') quantile(all.mse_theta, 0.025) quantile(all.mse_theta, 0.975)]; 
    % Mean simulated correlation between S and C
    res.corr_S_C = mean(all.corr_S_C, 'omitnan');
    % Mean simulated correlation between S and Y
    res.corr_S_Y = mean(all.corr_S_Y, 'omitnan');
    % Mean simulated correlation between C and Y
    res.corr_C_Y = mean(all.corr_C_Y, 'omitnan');

    disp(res);  % Display results

    % Save simulation summary results
    if length(samp_n) > length(iter)  % If aggregating over samples instead of iterations
        save(strcat(out_dir, 'summary_', model, '_scen', num2str(scenario), '_iter', num2str(iter), '_sample'), 'all', 'res');
    else                              % Aggregating over iterations
        save(strcat(out_dir, 'summary_', model, '_scen', num2str(scenario), '_samp', num2str(samp_n)), 'all', 'res');
    end    
end 