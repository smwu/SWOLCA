%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing code for weighted supervised OFMM and supervised OFMM       
% Programmer: SW             
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%54%%%%%%%%%%%%%%

% Test simulation code
% Population
load("simdata_scen2_iter2.mat");
tabulate(sim_data.true_Si);
tabulate(sim_data.true_Ci);  % true dist of pi
tabulate(sim_data.true_Ci(sim_data.true_Si==1));
tabulate(sim_data.Y_data(sim_data.true_Si==1 & sim_data.true_Ci==1));
tabulate(sim_data.Y_data(sim_data.true_Si==1 & sim_data.true_Ci==2));
sum(sim_data.X_data==2 & sim_data.true_Ci==1)
% Population thetas
sub_ci_1 = find(sim_data.true_Ci==1);
tabulate(sim_data.X_data(sub_ci_1,1));
tabulate(sim_data.X_data(sub_ci_1,2));
tabulate(sim_data.X_data(sub_ci_1,3));
tabulate(sim_data.X_data(sub_ci_1,4));
sub_ci_2 = find(sim_data.true_Ci==2);
tabulate(sim_data.X_data(sub_ci_2,1));
tabulate(sim_data.X_data(sub_ci_2,2));
tabulate(sim_data.X_data(sub_ci_2,3));
tabulate(sim_data.X_data(sub_ci_2,4));
corr(sim_data.true_Si, sim_data.true_Ci)
corr(sim_data.true_Si, sim_data.Y_data)
corr(sim_data.true_Ci, sim_data.Y_data)

% Sample
load("simdata_scen14_iter2_samp1.mat");
tabulate(sim_data.true_Si);
tabulate(sim_data.true_Ci);
tabulate(sim_data.true_Ci(sim_data.true_Si==1));
tabulate(sim_data.Y_data(sim_data.true_Si==1 & sim_data.true_Ci==1));
tabulate(sim_data.Y_data(sim_data.true_Si==1 & sim_data.true_Ci==2));
sum((sim_data.true_Si==1) .* sim_data.sample_wt)
sum((sim_data.true_Ci==1) .* sim_data.sample_wt)
sum(sim_data.X_data==2 & sim_data.true_Ci==1)
sum((sim_data.X_data==2 & sim_data.true_Ci==1) .* sim_data.sample_wt)
% Sample thetas
sub_ci_1 = find(sim_data.true_Ci==1);
tabulate(sim_data.X_data(sub_ci_1,1));
tabulate(sim_data.X_data(sub_ci_1,2));
tabulate(sim_data.X_data(sub_ci_1,3));
tabulate(sim_data.X_data(sub_ci_1,4));
sub_ci_2 = find(sim_data.true_Ci==2);
tabulate(sim_data.X_data(sub_ci_2,1));
tabulate(sim_data.X_data(sub_ci_2,2));
tabulate(sim_data.X_data(sub_ci_2,3));
tabulate(sim_data.X_data(sub_ci_2,4));
corr(sim_data.true_Si, sim_data.true_Ci)
corr(sim_data.true_Si, sim_data.Y_data)
corr(sim_data.true_Ci, sim_data.Y_data)

load('wsOFMM_latent_results_scen14_iter2_samp1.mat');
analysis.pi_med
tabulate(analysis.c_i)
analysis.theta_med
analysis.xi_med
normcdf(analysis.xi_med(1)+analysis.xi_med(3))
normcdf(analysis.xi_med(1)+analysis.xi_med(4))
normcdf(analysis.xi_med(2)+analysis.xi_med(3))
normcdf(analysis.xi_med(2)+analysis.xi_med(4))

load('sOFMM_latent_results_scen14_iter2_samp1.mat');
analysis.pi_med
tabulate(analysis.c_i)
analysis.theta_med
analysis.xi_med
normcdf(analysis.xi_med(1)+analysis.xi_med(3))
normcdf(analysis.xi_med(1)+analysis.xi_med(4))
normcdf(analysis.xi_med(2)+analysis.xi_med(3))
normcdf(analysis.xi_med(2)+analysis.xi_med(4))

% Time one run of sOFMM
tic
sOFMM_main(3,22);
toc
% timeit(@() sOFMM_main(3,22));

%%%%%%%% Loops
% Simulate population datasets for Scenarios 1,2 for 50 iterations 
for i = 1:50
    disp(i);
    sim_precision(i);
    sim_confounder(i);
end    

% Simulate samples for Scenario 2 iter 1, for 20 sampled sets each
iter = 1;
scen = 2;
for samp = 1:20
    sample_SRS(scen, iter, samp);
    sample_strat_eq(scen, iter, samp);
end    

% Simulate samples for Scenarios 1,2 iter 1-50, for 1 sampled set each
samp = 1;
for scen = 1:2
    for iter = 1:50
        sample_strat_eq(scen, iter, samp);
    end    
end

% Run sOFMM and wsOFMM models for simulated population data
tic
for scen = 1:2
    for iter = 1:50
        sOFMM_main(scen, iter, 0); 
        wsOFMM_main(scen, iter, 0); 
        wsOFMM_main_latent(scen, iter, 0); 
    end 
end    
toc

tic
% Run sOFMM model and wsOFMM model for simulated sample data
samp = 1;
for scen = 13:14
    for iter = 1:50
        sOFMM_main(scen, iter, samp);
        wsOFMM_main(scen, iter, samp);
        wsOFMM_latent_main(scen, iter, samp);
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


subset11 = find(sim_data.true_Si==1 & sim_data.true_Ci==1);
tabulate(sim_data.X_data(subset11))

subset1 = find(data_vars.wt_kappa==0.25);
tabulate(OFMM_params.c_i(subset1))