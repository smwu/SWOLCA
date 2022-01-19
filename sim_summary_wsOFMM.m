%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Summary Results    %
% for Weighted Supervised OFMM  %
% Programmer: SW                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;
% num_sims = 500;  % Number of simulated sets
num_sims = 1;
scen = 1;

% Initialize outputs, stored for all simulated sets
K_all = zeros(num_sims, 1);         % Initialize estimated number of classes
K_match_all = zeros(num_sims, 1);   % Initialize boolean vector of correct estimation of number of classes
sens_all = zeros(num_sims, 1);      % Initialize sensitivity
spec_all = zeros(num_sims, 1);      % Initialize specificity
dic6_all = zeros(num_sims, 1);      % Initialize DIC version 6 to penalize overfitting 
mse_Phi_all = zeros(num_sims, 1);   % Initialize MSE of class-specific response proportions
mse_theta_all = zeros(num_sims, 1); % Initialize MSE of class-specific item response probs

for sim_n = 1:num_sims              % For each simulated set
    load(strcat('simdata_wsRPC_scen', num2str(scen), '_iter', num2str(sim_n), '.mat'))  % Load simulated dataset
    load(strcat('wsOFMM_Results_scen', num2str(scen), '_iter', num2str(sim_n), '.mat')) % Load wsOFMM results
    
    K_all(sim_n) = analysis.k_red;                          % Estimated number of classes
    K_match_all(sim_n) = analysis.k_red == sim_data.true_K; % Estimated and true number of classes match (1/0)
    
    % Obtain sensitivity and specificity
    pw_classes_sens = zeros(analysis.k_red, sim_data.true_K); % Initialize matrix of pairwise class assigns for sensitivity
    pw_classes_spec = zeros(analysis.k_red, sim_data.true_K); % Initialize matrix of pairwise class assigns for specificity
    for j = 1:sim_data.true_K                                 % For each true class
        for i = 1:analysis.k_red                              % For each predicted class
            % Obtain prop subjects in same class who are correctly predicted to share a class
            pw_classes_sens(i, j) = sum(analysis.c_i == i & sim_data.true_Ci == j) / sum(sim_data.true_Ci == j);
            % Obtain prop subjects NOT in same class who are correctly predicted to NOT share a class
            pw_classes_spec(i, j) = sum(analysis.c_i ~= i & sim_data.true_Ci ~= j) / sum(sim_data.true_Ci == j);
        end
    end
    sens_class = max(pw_classes_sens);  % Max per column gives best pairing of true and predicted class membership    
    spec_class = max(pw_classes_spec);  % Max per column gives best pairing of true and predicted class membership 
    sens_all(sim_n) = mean(sens_class); % Mean sensitivity over all true classes
    spec_all(sim_n) = mean(spec_class); % Mean specificity over all true classes
    
    dic6_all(sim_n) = analysis.dic6;  % DIC version 6 to penalize overfitting
        
    % MSE of class-specific response proportions
    mse_Phi_all(sim_n) = immse(analysis.Phi_med, sim_data.true_Phi);  
    
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
    mse_theta_all(sim_n) = mean(mse_theta_class); % Total MSE is mean over all classes    
end

% Median and 95% CI of estimated number of classes over all simulated datasets
res.K_median = [median(K_all) quantile(K_all, 0.025, 2) quantile(K_all, 0.975, 2)];          
% Proportion of simulations with number of classes correctly estimated
res.K_match = mean(K_match_all);                     
% Median and 95% CI of sensitivity over all simulated datasets
res.sens_median = [median(sens_all) quantile(sens_all, 0.025, 2) quantile(sens_all, 0.975, 2)];  
% Median and 95% CI of specificity over all simulated datasets
res.spec_median = [median(spec_all) quantile(spec_all, 0.025, 2) quantile(spec_all, 0.975, 2)];
% Median and 95% CI of DIC-6 over all simulated datasets
res.dic6_median = [median(dic6_all, 'omitnan') quantile(dic6_all, 0.025, 2) quantile(dic6_all, 0.975, 2)];      
% Median and 95% CI of MSE of class-specific response proportions over all simulated datasets
res.mse_Phi_median = [median(mse_Phi_all, 'omitnan') quantile(mse_Phi_all, 0.025, 2) quantile(mse_Phi_all, 0.975, 2)]; 
% Median and 95% CI of MSE of item response probabilities over all simulated datasets 
res.mse_theta_median = [median(mse_theta_all) quantile(mse_theta_all, 0.025, 2) quantile(mse_theta_all, 0.975, 2)];       

% Save simulation summary results
save(strcat('sim_summary_scen', num2str(scen)), 'res');

%% Testing code
% num_sims = 1;
