%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing code for weighted supervised OFMM and supervised OFMM       
% Programmer: SW             
% 
% Assume equal subpopulation sizes
% Scenario 1: Full population. All weights equal to 1  
% Scenario 2: Sample 5% of total population (SRS). All weights equal
% Scenario 3: Sample 5% from each subpop (proportional allocation). 
%             All weights equal up to rounding      
% Scenario 4: Sample 1000 from each subpop (equal allocation). 
%             Diff weights per subpop 
%
% Data description:
% We assume individuals come from 4 subpopulations. Sampling of 4000 
% subjects is stratified by subpop. There are 3 global dietary patterns, 
% and for each subpopulation, there are 2 local dietary patterns. 
% Each food item for an individual belongs to a global pattern or to a 
% subpopulation-specific local pattern.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%54%%%%%%%%%%%%%%

%% Testing code
%%%%%%% For one iteration
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

% Check correlations
disp(corr(sim_data.true_Si, sim_data.true_Ci));
disp(corr(sim_data.true_Si, sim_data.Y_data));
disp(corr(sim_data.true_Ci, sim_data.Y_data)); 

% Time one run of sOFMM
tic
sOFMM_main(3,22);
toc
% timeit(@() sOFMM_main(3,22));

%%%%%%%% Loops
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
    
% Display summaries for multiple population scenarios
for scen = 1:4
    sim_summary_wsOFMM(scen, 10, "wsOFMM")
end

% Display summaries for multiple sample scenarios
for scen = 5:16
    sim_summary_wsOFMM_sample(scen, 1:50, 1, "wsOFMM")
end


% Display correlations for multiple simulated population datasets and multiple
% scenarios
num_iters = 50;
data_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Data/"; 
for scen = 1:4
    all.corr_S_C = NaN(num_iters, 1);  % Initialize correlation between S and C
    all.corr_S_Y = NaN(num_iters, 1);  % Initialize correlation between S and Y
    all.corr_C_Y = NaN(num_iters, 1);  % Initialize correlation between C and Y
    for sim_n = 1:num_iters % For each simulated set or sample

        % Load simulated dataset
        load(strcat(data_dir, 'simdata_scen', num2str(scen), '_iter', num2str(sim_n), '.mat'), 'sim_data') 

        % Get simulated correlations
        all.corr_S_C(sim_n) = corr(sim_data.true_Si, sim_data.true_Ci);
        all.corr_S_Y(sim_n) = corr(sim_data.true_Si, sim_data.Y_data);
        all.corr_C_Y(sim_n) = corr(sim_data.true_Ci, sim_data.Y_data);

    end
    disp(strcat('Scenario ', num2str(scen)));
    disp(mean(all.corr_S_C));
    disp(mean(all.corr_S_Y));
    disp(mean(all.corr_C_Y));
    disp(strcat('sd S C', num2str(std(all.corr_S_C))));
end    

% Display correlations for multiple simulated sample datasets and multiple
% scenarios
num_iters = 50;
samp_iter = 1;
data_dir = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/Data/"; 
for scen = 5:16
    all.corr_S_C = NaN(num_iters, 1);  % Initialize correlation between S and C
    all.corr_S_Y = NaN(num_iters, 1);  % Initialize correlation between S and Y
    all.corr_C_Y = NaN(num_iters, 1);  % Initialize correlation between C and Y
    for sim_n = 1:num_iters % For each simulated set or sample

        % Load simulated dataset
        load(strcat(data_dir, 'simdata_scen', num2str(scen), '_iter', num2str(sim_n), '_samp', num2str(samp_iter), '.mat'), 'sim_data') 

        % Get simulated correlations
        all.corr_S_C(sim_n) = corr(sim_data.true_Si, sim_data.true_Ci);
        all.corr_S_Y(sim_n) = corr(sim_data.true_Si, sim_data.Y_data);
        all.corr_C_Y(sim_n) = corr(sim_data.true_Ci, sim_data.Y_data);

    end
    disp(strcat('Scenario ', num2str(scen)));
    disp(mean(all.corr_S_C));
    disp(mean(all.corr_S_Y));
    disp(mean(all.corr_C_Y));
    disp(strcat('sd S C', num2str(std(all.corr_S_C))));
end   