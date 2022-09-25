%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Summary Results for wOFMM bootPRS
% Programmer: SW                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_reps = 100;  % Number of simulated samples
n_PRS = 200;   % Number of PRS's
%n_reps = 2;
%n_PRS = 2;
k_max = 50;
res_dir = strcat(pwd, "/M200/"); % Directory with model results
true_dir = strcat(pwd, "/");  % Where the true simulated population is
%res_dir = pwd;   

% Initialize outputs, stored for all simulated sets
all_mse_pi_mean = NaN(n_reps, 1);   % Initialize MSE of pi with mean
all_mse_pi_med = NaN(n_reps, 1);    % Initialize MSE of pi with median
all_mse_theta0 = NaN(n_reps, 1); % Initialize MSE of theta0
all_runtime = NaN(n_reps, 1); % Initialize runtime
pi_keepwt = NaN(n_reps, 3); % Initialize posterior pi for boxplots, assuming k_true=3

for sim_n = 1:n_reps      % For each simulated set or sample
    
    disp(strcat("sim_n: ",num2str(sim_n)));
    
    % Initialize pooled parameter vectors
    pi_mean_PRS = NaN(n_PRS, k_max);
    pi_med_PRS = NaN(n_PRS, k_max);

    % Obtain posterior estimate pooling across PRS samples
    for PRS_iter = 1:n_PRS
        if mod(PRS_iter, 50) == 0
            disp(strcat("PRS_iter: ", num2str(PRS_iter)));
        end
        % Construct file names                             
        file_str = strcat(res_dir, 'bwsOFMM_latent_results_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n), '_prs', num2str(PRS_iter), '.mat');

        % Check if the results file exists
        if isfile(file_str) 
            % Load results
            load(file_str);

            % Mean parameter estimates across MCMC samples
            theta0_mean = reshape(mean(theta0_order), [p,k0,d_max]);
            pi_mean = mean(pi_order);
            [val0, ind0] = max(theta0_mean, [], 3); % Modal pattern of global clusters
            [t0, ia, ic] = unique(ind0', 'rows');
            k_red_mean = length(ia);  % Number of unique clusters
            pi_mean_PRS(PRS_iter, 1:k_red_mean) = pi_mean(ia) / sum(pi_mean(ia));  
            theta0_mean = theta0_mean(:, ia, :);
            [maxv_theta_mean, maxp_theta_mean] = max(theta0_mean, [], 3);

            % Median parameter estimates across MCMC samples
            theta0_med = reshape(median(theta0_order),[p,k0,d_max]);
            pi_med = median(pi_order);
            [val0, ind0] = max(theta0_med, [], 3); % Modal patterns
            [t0, ia, ic] = unique(ind0', 'rows');
            k_red_med = length(ia);  % Number of unique clusters
            pi_med_PRS(PRS_iter, 1:k_red_med) = pi_med(ia) / sum(pi_med(ia));
            theta0_med = theta0_med(:, ia, :);
            [maxv_theta_med, maxp_theta_med] = max(theta0_med, [], 3);

        else
            % If the results file does not exist, print out statement
            disp(strcat('Iteration ', num2str(sim_n), ' population ', num2str(PRS_iter), ' is missing.'));
        end    
    end
    
    % Load true data
    load(strcat(true_dir, '/uOFMsimdata_A',num2str(sim_n),'.mat'))

    % Get MSE for pi combining MI runs
    pi_mean_all = mean(pi_mean_PRS, 1, 'omitnan');  % mean of mean pi's
    pi_mean_all = pi_mean_all(~isnan(pi_mean_all)); % Remove NaN's
    pi_mean_all = pi_mean_all / sum(pi_mean_all);   % Normalize
    
    pi_med_all = mean(pi_med_PRS, 1, 'omitnan'); 
    pi_med_all = pi_med_all(~isnan(pi_med_all));
    pi_med_all = pi_med_all / sum(pi_med_all);
    
    %     k_true = length(true_pi); 
    %     k_pred = length(pi_mean_all);
    %     test_true=zeros(k_pred,k_true);
    %     for ti = 1:k_true
    %         for zi = 1:k_pred
    %             test_true(zi, ti)=sum(z_max == zi & sample_Ci == ti);
    %         end
    %         
    %     end
    %     [v_ind,t_ind]=max(test_true);
    %     tablCi_sample=tabulate(sample_Ci);
    %     classwt_err(iset,:)=v_ind./transpose(tablCi_sample(:,2));
    %     pi_keepwt(iset,:)=pi_med(t_ind);
    %     pi_mse=immse(pi_med(t_ind),true_pi);
    
    % Get optimal order of true_pi to match pi_posterior
    % Assumes number of clusters matches
    if length(pi_mean_all) == length(true_pi)
        perms_K = perms(1:3);          % permutations of (1,2,3). Total 3! combos
        mse_mean_mat = zeros(factorial(3),1);
        mse_med_mat = zeros(factorial(3), 1);
        for pp = 1:factorial(3)
            pi_perm = true_pi([perms_K(pp,:)]);
            mse_mean_mat(pp) = immse(pi_mean_all, pi_perm);
            mse_med_mat(pp) = immse(pi_med_all, pi_perm);
        end
        [pi_mean_mse, ind_mean] = min(mse_mean_mat);
        [pi_med_mse, ind_med] = min(mse_med_mat);
        pi_keepwt(sim_n, :) = pi_mean_all;
    else  % When true pi and estimated pi have different number of components
        N_true = length(true_pi); 
        N_est = length(pi_mean_all);
        if N_est > N_true  % if estimated pi has more elements than true pi
            pi_red = pi_mean_all;
            pi_fixed = true_pi;
            N = N_est;
            K = N_true;
        else
            pi_red = true_pi;
            pi_fixed = pi_mean_all;
            N = N_true;
            K = N_est;
        end    
        % Get all permutations selection K out of N elements
        perms_K = nchoosek(1:N,K);
        perms_K = reshape(perms_K(:,perms(1:K)),[],K);
        n_perm = length(perms_K);
        % Get optimal order of true_pi to match pi_est
        % when number of clusters differs
        mse_mean_mat = zeros(n_perm,1);
        mse_med_mat = zeros(n_perm, 1);
        for pp = 1:n_perm
            pi_perm = pi_red([perms_K(pp,:)]);
            mse_mean_mat(pp) = immse(pi_fixed, pi_perm);
            mse_med_mat(pp) = immse(pi_fixed, pi_perm);
        end
        [pi_mean_mse, ind_mean] = min(mse_mean_mat);
        [pi_med_mse, ind_med] = min(mse_med_mat);
        if N_est > N_true
            pi_keepwt(sim_n, :) = pi_mean_all(perms_K(ind_mean,:));
        else
            pi_keepwt(sim_n, :) = pi_mean_all;
        end    

    end    
    %pi_perm_mean = true_pi([perms_K(ind_mean, :)]);
    %pi_perm_med = true_pi([perms_K(ind_med, :)]);
    
    all_mse_pi_mean(sim_n) = pi_mean_mse;
    all_mse_pi_med(sim_n) = pi_med_mse;
    all_runtime(sim_n) = eltime;
end

% Mean and 95% CI over all iterations 
mean_mse_pi = [mean(all_mse_pi_mean, 'omitnan') quantile(all_mse_pi_mean, 0.025) quantile(all_mse_pi_mean, 0.975)];          
med_mse_pi = [mean(all_mse_pi_med, 'omitnan') quantile(all_mse_pi_med, 0.025) quantile(all_mse_pi_med, 0.975)];          
mean_runtime = [mean(all_runtime, 'omitnan') quantile(all_runtime, 0.025) quantile(all_runtime, 0.975)];

disp("pi (mean)");
disp(mean_mse_pi);  % Display results
disp("pi (median)");
disp(med_mse_pi);
disp("runtime");
disp(mean_runtime);

save(strcat(res_dir, 'res_bootPRS_A', num2str(sim_n)), 'mean_mse_pi', 'mean_runtime', 'all_mse_pi_mean', 'pi_keepwt');

figure;
boxplot(sort(pi_keepwt,2))
ylim([0 0.7])
yline(0.6,'LineStyle','-.','Color','magenta','LineWidth',1)
yline(0.3,'LineStyle','-.','Color','green','LineWidth',1)
yline(0.1,'LineStyle','-.','Color','black','LineWidth',1)
saveas(gcf, 'bootPRS_A.png')

