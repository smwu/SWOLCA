%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Summary Results    %
% for Weighted Supervised OFMM  %
% Programmer: SW                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;
% num_sims = 500;  % Number of simulated sets
num_sims = 1;

% Initialize outputs, stored for all simulated sets
K_all = zeros(num_sims, 1);         % Initialize estimated number of classes
sens_all = zeros(num_sims, 1);      % Initialize sensitivity
dic_all = zeros(num_sims, 1);       % Initialize DIC 
dic6_all = zeros(num_sims, 1);      % Initialize DIC version 6 to penalize overfitting 
mse_Phi_all = zeros(num_sims, 1);       % Initialize MSE of class-specific response proportions
    % mse_theta_all = zeros(num_sims, 1); % Initialize MSE of class-specific item response probs
    % pi_all = zeros(num_sims, 1);        % Initialize class membership probabilities
    % theta_all = zeros(num_sims, 1);     % Initialize item response probabilities
    % xi_all = zeros(num_sims, 1);        % Initialize probit model coefficients
    % spec_all;

for iset = 1:num_sims               % For each simulated set
    load(strcat('simdata_wsRPC_scen4_iter', num2str(iset), '.mat'))  % Load simulated dataset
    %load(strcat('wsOFMM_Results_1_iter', num2str(iset), '.mat')) % Load wsOFMM results
    load(strcat('wsOFMM_Results_4_', num2str(iset), '.mat')) % Load wsOFMM results
    true_K = 3;

    K_all(iset) = analysis.k_red;  % Estimated number of classes
    
    % Obtain sensitivity
    pairwise_classes = zeros(analysis.k_red, true_K); % Initialize matrix of pairwise class assignments
    for j = 1:true_K                                  % for each true class
        for i = 1:analysis.k_red                      % for each predicted class
            % Obtain the proportion of subjects in the same class who are correctly predicted to share a class
            pairwise_classes(i, j) = sum(analysis.c_i == i & true_Ci == j) / sum(true_Ci == j);
        end
    end
    sens_class = max(pairwise_classes); % Max of each column, i.e., best pairing of true and predicted class membership       
    sens_all(iset) = mean(sens_class);  % Mean sensitivity over all true classes
    
    dic_all(iset) = analysis.dic;    % DIC
    dic6_all(iset) = analysis.dic; 
    % dic6_all(iset) = analysis.dic6;  % DIC version 6 to penalize overfitting
        
    % MSE of class-specific response proportions
    mse_Phi_all(iset) = immse(analysis.Phi_med, true_Phi);  
    %     % MSE of class-specific item response probabilities
    %     mse_theta_all(iset) = immse(analysis.theta_med, global_thetas); 
    
    %     pi_all(iset) = analysis.pi_med;        % Posterior median for cluster membership probs
    %     theta_all(iset) = analysis.theta_med;  % Posterior median for item response probs
    %     xi_all(iset) = analysis.xi_med;        % Posterior median for probit model coefs
    
end

% Print out median estimates over all simulated datasets
K_median = median(K_all)                        % Median estimated number of classes
sens_median = median(sens_all)                  % Median sensitivity
dic_median = median(dic_all, 'omitnan')         % Median DIC
dic6_median = median(dic6_all, 'omitnan')       % Median DIC-6
mse_Phi_median = median(mse_Phi_all, 'omitnan') % Median MSE of class-specific response proportions
    % mse_theta_median = median(mse_theta_all)        % Median MSE of item response probabilities
    % pi_median = median(pi_all)                      % Median class membership probabilities
    % theta_median = median(theta_all)                % Median item response probabilities
    % xi_median = median(xi_all)                      % Median probit model coefficients
