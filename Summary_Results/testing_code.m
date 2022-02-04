%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing code for weighted supervised OFMM and supervised OFMM       
% Programmer: SW             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Testing code
% for one iteration
in_dir = "/n/home01/stephwu18/wsOFMM/data/";   % Input directory
out_dir = "/n/home01/stephwu18/wsOFMM/data/";  % Output directory
scenario = 5;                                  % Scenario for full population
sim_n = 1;                                     % Simulation iteration
% Load population simulated dataset
load(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '.mat'), 'sim_data')  

% MSE for Phi
disp(immse(analysis.Phi_med, sim_data.true_Phi))
% Obtain sensitivity and specificity for sOFMM
sens = zeros(analysis.k_red, sim_data.true_K); % Initialize matrix of pairwise class assigns for sensitivity
spec = zeros(analysis.k_red, sim_data.true_K); % Initialize matrix of pairwise class assigns for specificity
for j = 1:sim_data.true_K                                 % For each true class
    for i = 1:analysis.k_red                              % For each predicted class
        % Obtain prop subjects in same class who are correctly predicted to share a class
        sens(i, j) = sum(analysis.c_i == i & sim_data.true_Ci == j) / sum(sim_data.true_Ci == j);
        % Obtain prop subjects NOT in same class who are correctly predicted to NOT share a class
        spec(i, j) = sum(analysis.c_i ~= i & sim_data.true_Ci ~= j) / sum(sim_data.true_Ci ~= j);
    end
end
sens_class = max(sens);  % Max per column gives best pairing of true and predicted class membership    
spec_class = max(spec);  % Max per column gives best pairing of true and predicted class membership 
disp(mean(sens_class)); % Mean sensitivity over all true classes
disp(mean(spec_class)); % Mean specificity over all true classes

% MSE of class-specific item response probabilities for sOFMM
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
disp(mean(mse_theta_class)); % Total MSE is mean over all classes

% Obtain sensitivity and specificity for sLCA
k_red = length(pi_med);
sens = zeros(k_red, sim_data.true_K); % Initialize matrix of pairwise class assigns for sensitivity
spec = zeros(k_red, sim_data.true_K); % Initialize matrix of pairwise class assigns for specificity
for j = 1:sim_data.true_K                                 % For each true class
    for i = 1:k_red                              % For each predicted class
        % Obtain prop subjects in same class who are correctly predicted to share a class
        sens(i, j) = sum(pred_ci == i & sim_data.true_Ci == j) / sum(sim_data.true_Ci == j);
        % Obtain prop subjects NOT in same class who are correctly predicted to NOT share a class
        spec(i, j) = sum(pred_ci ~= i & sim_data.true_Ci ~= j) / sum(sim_data.true_Ci ~= j);
    end
end
sens_class = max(sens);  % Max per column gives best pairing of true and predicted class membership    
spec_class = max(spec);  % Max per column gives best pairing of true and predicted class membership 
disp(mean(sens_class)); % Mean sensitivity over all true classes
disp(mean(spec_class)); % Mean specificity over all true classes

% MSE of class-specific item response probabilities for sLCA
pw_classes_mse = zeros(k_red, sim_data.true_K);         % Initialize matrix of pairwise class MSEs
dim = size(sim_data.true_global_thetas); p = dim(1); d = dim(3); % Get food data dimensions
for j = 1:sim_data.true_K                                        % For each true class
    % Obtain the true item response probs for true class j
    global_thetas_j = reshape(sim_data.true_global_thetas(:, j, :), [p, d]); 
    for i = 1:k_red                                     % For each predicted class
        % Obtain the predicted item response probs for predicted class i
        theta_med_i = reshape(theta0_med(:, i, :), [p, d]);
        % Obtain MSE of true and predicted item response probs for each pairwise class combo
        pw_classes_mse(i, j) = immse(theta_med_i, global_thetas_j);
    end
end
mse_theta_class = min(pw_classes_mse);        % Min per column gives best pairing of true and predicted class membership    
disp(mean(mse_theta_class)); % Total MSE is mean over all classes

% Time one run of sOFMM
tic
sOFMM_main(3,22);
toc
% timeit(@() sOFMM_main(3,22));

% Simulate datasets for Scenarios 1,2,3,4 for 10 iterations
for i = 1:10
    sim_no_local(i);
    sim_no_local_pi_subpop(i);
    sim_uneq(i);
    sim_uneq_pi_subpop(i);
end    

% Simulate samples for Scenarios 1,2,3,4 iter 1, for 5 sampled sets each
iter = 1;
for scen = 1:4
    for samp = 1:5
        sample_SRS(scen, iter, samp);
        sample_strat_prop(scen, iter, samp);
        sample_strat_eq(scen, iter, samp);
    end    
end

% Run sOFMM model for simulated population data
tic
for scen = 1:4
    for iter = 1:10
        sOFMM_main(1, iter, 0); 
    end 
end    
toc

tic
% Run sOFMM model and wsOFMM model for simulated sample data
iter = 1;
for scen = 1:4
    for samp = 1:5
        sOFMM_main(1, iter, samp);
        wsOFMM_main(1, iter, samp);
    end    
end
toc
    
% Check correlations
disp(corr(sim_data.true_Si, sim_data.true_Ci));
disp(corr(sim_data.true_Si, sim_data.Y_data));
disp(corr(sim_data.true_Ci, sim_data.Y_data)); 

