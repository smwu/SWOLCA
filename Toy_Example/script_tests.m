% Input directory
in_dir = "C:/Users/Lang/Documents/Harvard/Research/Briana/supRPC/wsOFMM/Toy_Example/";   
% Output directory 
out_dir = "C:/Users/Lang/Documents/Harvard/Research/Briana/supRPC/wsOFMM/Toy_Example/";  
% Population scenario 2 (confounder), iteration 1
pop_data = importdata(strcat(in_dir, 'simdata_scen2_iter1.mat'));
% Sample scenario 14 (unequal sampling), iteration 1, sample 1
samp_data = importdata(strcat(in_dir, 'simdata_scen14_iter1_samp1.mat'));

% Define an absolute tolerance
tol = 1e-4;

% Note: tests to be run SEQUENTIALLY. Later tests depend on the output of 
% earlier tests

%% Observed population characteristics
% P(S_i=1)
sum(pop_data.true_Si==1) / numel(pop_data.true_Si)
% P(C_i=1)
sum(pop_data.true_Ci==1) / numel(pop_data.true_Ci)
% P(C_i=1 | S_i=1)
sum(pop_data.true_Ci==1 & pop_data.true_Si==1) / sum(pop_data.true_Si==1)
% P(Y_i=1 | S_i=1, C_i=1)
sum(pop_data.Y_data==1 & pop_data.true_Ci==1 & pop_data.true_Si==1) / ...
    sum(pop_data.true_Ci==1 & pop_data.true_Si==1)
% P(Y_i=1 | S_i=1, C_i=2)
sum(pop_data.Y_data==1 & pop_data.true_Ci==2 & pop_data.true_Si==1) / ...
    sum(pop_data.true_Ci==2 & pop_data.true_Si==1)
% P(Y_i=1 | S_i=2, C_i=1)
sum(pop_data.Y_data==1 & pop_data.true_Ci==1 & pop_data.true_Si==2) / ...
    sum(pop_data.true_Ci==1 & pop_data.true_Si==2)
% P(Y_i=1 | S_i=2, C_i=2)
sum(pop_data.Y_data==1 & pop_data.true_Ci==2 & pop_data.true_Si==2) / ...
    sum(pop_data.true_Ci==2 & pop_data.true_Si==2)
% P(X_i1=1 | C_i=1)
sum(pop_data.X_data(:,1)==1 & pop_data.true_Ci==1) / sum(pop_data.true_Ci==1)
% P(X_i2=1 | C_i=1)
sum(pop_data.X_data(:,2)==1 & pop_data.true_Ci==1) / sum(pop_data.true_Ci==1)
% P(X_i3=1 | C_i=1)
sum(pop_data.X_data(:,3)==1 & pop_data.true_Ci==1) / sum(pop_data.true_Ci==1)
% P(X_i4=1 | C_i=1)
sum(pop_data.X_data(:,4)==1 & pop_data.true_Ci==1) / sum(pop_data.true_Ci==1)
% P(X_i1=1 | C_i=2)
sum(pop_data.X_data(:,1)==1 & pop_data.true_Ci==2) / sum(pop_data.true_Ci==2)
% P(X_i2=1 | C_i=2)
sum(pop_data.X_data(:,2)==1 & pop_data.true_Ci==2) / sum(pop_data.true_Ci==2)
% P(X_i3=1 | C_i=2)
sum(pop_data.X_data(:,3)==1 & pop_data.true_Ci==2) / sum(pop_data.true_Ci==2)
% P(X_i4=1 | C_i=2)
sum(pop_data.X_data(:,4)==1 & pop_data.true_Ci==2) / sum(pop_data.true_Ci==2)
% Cor(S,C)
corr(pop_data.true_Si, pop_data.true_Ci)
% Cor(C,Y)
corr(pop_data.true_Ci, pop_data.Y_data)
% Cor(S,Y)
corr(pop_data.true_Si, pop_data.Y_data)
% P(X_i1=1 | S_i=1, C_i=1)
sum(pop_data.X_data(:,1)==1 & pop_data.true_Ci==1 & pop_data.true_Si==1) / ...
    sum(pop_data.true_Ci==1 & pop_data.true_Si==1)
% P(X_i1=1 | S_i=2, C_i=1)
sum(pop_data.X_data(:,1)==1 & pop_data.true_Ci==1 & pop_data.true_Si==2) / ...
    sum(pop_data.true_Ci==1 & pop_data.true_Si==2)
% P(X_i1=1 | S_i=1, C_i=2)
sum(pop_data.X_data(:,1)==1 & pop_data.true_Ci==2 & pop_data.true_Si==1) / ...
    sum(pop_data.true_Ci==2 & pop_data.true_Si==1)
% P(X_i1=1 | S_i=2, C_i=2)
sum(pop_data.X_data(:,1)==1 & pop_data.true_Ci==2 & pop_data.true_Si==2) / ...
    sum(pop_data.true_Ci==2 & pop_data.true_Si==2)
% P(X_i1=2 | S_i=1, C_i=1)
sum(pop_data.X_data(:,2)==1 & pop_data.true_Ci==1 & pop_data.true_Si==1) / ...
    sum(pop_data.true_Ci==1 & pop_data.true_Si==1)
% P(X_i1=2 | S_i=2, C_i=1)
sum(pop_data.X_data(:,2)==1 & pop_data.true_Ci==1 & pop_data.true_Si==2) / ...
    sum(pop_data.true_Ci==1 & pop_data.true_Si==2)
% P(X_i1=2 | S_i=1, C_i=2)
sum(pop_data.X_data(:,2)==1 & pop_data.true_Ci==2 & pop_data.true_Si==1) / ...
    sum(pop_data.true_Ci==2 & pop_data.true_Si==1)
% P(X_i1=2 | S_i=2, C_i=2)
sum(pop_data.X_data(:,2)==1 & pop_data.true_Ci==2 & pop_data.true_Si==2) / ...
    sum(pop_data.true_Ci==2 & pop_data.true_Si==2)

%% Observed sample characteristics
% P(S_i=1)
sum(samp_data.true_Si==1) / numel(samp_data.true_Si)
% P(C_i=1)
sum(samp_data.true_Ci==1) / numel(samp_data.true_Ci)
% P(C_i=1 | S_i=1)
sum(samp_data.true_Ci==1 & samp_data.true_Si==1) / sum(samp_data.true_Si==1)
% P(Y_i=1 | S_i=1, C_i=1)
sum(samp_data.Y_data==1 & samp_data.true_Ci==1 & samp_data.true_Si==1) / ...
    sum(samp_data.true_Ci==1 & samp_data.true_Si==1)
% P(Y_i=1 | S_i=1, C_i=2)
sum(samp_data.Y_data==1 & samp_data.true_Ci==2 & samp_data.true_Si==1) / ...
    sum(samp_data.true_Ci==2 & samp_data.true_Si==1)
% P(Y_i=1 | S_i=2, C_i=1)
sum(samp_data.Y_data==1 & samp_data.true_Ci==1 & samp_data.true_Si==2) / ...
    sum(samp_data.true_Ci==1 & samp_data.true_Si==2)
% P(Y_i=1 | S_i=2, C_i=2)
sum(samp_data.Y_data==1 & samp_data.true_Ci==2 & samp_data.true_Si==2) / ...
    sum(samp_data.true_Ci==2 & samp_data.true_Si==2)
% P(X_i1=1 | C_i=1)
sum(samp_data.X_data(:,1)==1 & samp_data.true_Ci==1) / sum(samp_data.true_Ci==1)
% P(X_i2=1 | C_i=1)
sum(samp_data.X_data(:,2)==1 & samp_data.true_Ci==1) / sum(samp_data.true_Ci==1)
% P(X_i3=1 | C_i=1)
sum(samp_data.X_data(:,3)==1 & samp_data.true_Ci==1) / sum(samp_data.true_Ci==1)
% P(X_i4=1 | C_i=1)
sum(samp_data.X_data(:,4)==1 & samp_data.true_Ci==1) / sum(samp_data.true_Ci==1)
% P(X_i1=1 | C_i=2)
sum(samp_data.X_data(:,1)==1 & samp_data.true_Ci==2) / sum(samp_data.true_Ci==2)
% P(X_i2=1 | C_i=2)
sum(samp_data.X_data(:,2)==1 & samp_data.true_Ci==2) / sum(samp_data.true_Ci==2)
% P(X_i3=1 | C_i=2)
sum(samp_data.X_data(:,3)==1 & samp_data.true_Ci==2) / sum(samp_data.true_Ci==2)
% P(X_i4=1 | C_i=2)
sum(samp_data.X_data(:,4)==1 & samp_data.true_Ci==2) / sum(samp_data.true_Ci==2)
% Cor(S,C)
corr(samp_data.true_Si, samp_data.true_Ci)
% Cor(C,Y)
corr(samp_data.true_Ci, samp_data.Y_data)
% Cor(S,Y)
corr(samp_data.true_Si, samp_data.Y_data)
% P(X_i1=1 | S_i=1, C_i=1)
sum(samp_data.X_data(:,1)==1 & samp_data.true_Ci==1 & samp_data.true_Si==1) / ...
    sum(samp_data.true_Ci==1 & samp_data.true_Si==1)
% P(X_i1=1 | S_i=2, C_i=1)
sum(samp_data.X_data(:,1)==1 & samp_data.true_Ci==1 & samp_data.true_Si==2) / ...
    sum(samp_data.true_Ci==1 & samp_data.true_Si==2)
% P(X_i1=1 | S_i=1, C_i=2)
sum(samp_data.X_data(:,1)==1 & samp_data.true_Ci==2 & samp_data.true_Si==1) / ...
    sum(samp_data.true_Ci==2 & samp_data.true_Si==1)
% P(X_i1=1 | S_i=2, C_i=2)
sum(samp_data.X_data(:,1)==1 & samp_data.true_Ci==2 & samp_data.true_Si==2) / ...
    sum(samp_data.true_Ci==2 & samp_data.true_Si==2)
% P(X_i1=2 | S_i=1, C_i=1)
sum(samp_data.X_data(:,2)==1 & samp_data.true_Ci==1 & samp_data.true_Si==1) / ...
    sum(samp_data.true_Ci==1 & samp_data.true_Si==1)
% P(X_i1=2 | S_i=2, C_i=1)
sum(samp_data.X_data(:,2)==1 & samp_data.true_Ci==1 & samp_data.true_Si==2) / ...
    sum(samp_data.true_Ci==1 & samp_data.true_Si==2)
% P(X_i1=2 | S_i=1, C_i=2)
sum(samp_data.X_data(:,2)==1 & samp_data.true_Ci==2 & samp_data.true_Si==1) / ...
    sum(samp_data.true_Ci==2 & samp_data.true_Si==1)
% P(X_i1=2 | S_i=2, C_i=2)
sum(samp_data.X_data(:,2)==1 & samp_data.true_Ci==2 & samp_data.true_Si==2) / ...
    sum(samp_data.true_Ci==2 & samp_data.true_Si==2)
% Weighted P(S_i=1)
sum((samp_data.true_Si==1) .* samp_data.sample_wt) / numel(pop_data.true_Si)
% Weighted P(C_i=1)
sum((samp_data.true_Ci==1) .* samp_data.sample_wt) / numel(pop_data.true_Ci)

%% Test 1: wtd_get_data_vars function
data_vars = wtd_get_data_vars(samp_data);
assert(data_vars.d_max == 2)
assert(sum(data_vars.wt_kappa) == 200)
assert(data_vars.n == 200)
assert(data_vars.p == 4)
assert(isequal(data_vars.lin_idx(1:5), [1; 1; 5; 1; 1]))
assert(isequal(data_vars.wt_kappa_mat(1:2,1:2), [0.25 0.25; 1.75 1.75]))
assert(sum(data_vars.food, 'all') == 1180)
assert(sum(data_vars.y) == 104)

%% Test 2: wtd_init_OFMM_params function
k_max = 50;
sp_k = k_max;
alpha = ones(1, k_max) / sp_k;
eta = ones(1, data_vars.d_max);
rng(1, 'twister');
OFMM_params = wtd_init_OFMM_params(data_vars, k_max, alpha, eta);
assert(abs(OFMM_params.pi(7) - 0.2319) <= tol)
assert(abs(OFMM_params.theta(1,1,1) - 0.2654) <= tol)
assert(abs(OFMM_params.theta(4,50,2) - 0.0082) <= tol)
assert(isequal(OFMM_params.c_i(1:5), [48; 8; 29; 35; 48]))
assert(abs(OFMM_params.n_ci(7) - 37.5) <= tol)

%% Test 3: wtd_init_probit_params function
q_dem = dummyvar(samp_data.true_Si);        % Matrix of demographic covariates in cell-means format. Default contains subpop.                               
S = size(q_dem, 2);                         % Number of demographic covariates in the probit model
p_cov = k_max + S;                          % Number of covariates in probit model
mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVGamma(shape=5/2, scale=5/2)
Sig_0 = diag(Sig_0); 
probit_params = init_probit_params(data_vars, k_max, q_dem, mu_0, Sig_0);
assert(abs(probit_params.xi(1) - 1.7905) <= tol)
assert(abs(probit_params.xi(52) + 0.1900) <= tol)
assert(abs(probit_params.indiv_lik_probit_class(1,50) - 0.0547) <= tol)
assert(abs(probit_params.indiv_lik_probit_class(200,1) - 1) <= tol)

%% Test 4: wtd_update_MCMC function
n_runs = 9;
burn = 0;
thin = 1;
MCMC_out.pi = zeros(n_runs / thin, k_max);  
MCMC_out.theta = zeros(n_runs / thin, data_vars.p, k_max, data_vars.d_max); 
MCMC_out.c_i = zeros(n_runs / thin, data_vars.n);  
MCMC_out.xi = zeros(n_runs / thin, p_cov);  
MCMC_out.z_i = zeros(n_runs / thin, data_vars.n);
MCMC_out.loglik = zeros(n_runs / thin, 1); 
for iter = 1:n_runs
    [MCMC_out, OFMM_params, probit_params] = wtd_update_MCMC(MCMC_out, data_vars, ...
        OFMM_params, probit_params, thin, k_max, q_dem, alpha, eta, mu_0, Sig_0, iter);
end  
% n_runs=1;
% assert(abs(OFMM_params.pi(35) - 0.1315) <= tol)
% assert(abs(OFMM_params.theta(1,1,1) - 0.6526) <= tol)
% assert(abs(OFMM_params.theta(4,50,2) - 0.8460) <= tol) 
% assert(isequal(OFMM_params.c_i(1:5), [22; 24; 48; 7; 7]))
% assert(OFMM_params.n_ci(7) == 37.5)
% assert(abs(probit_params.xi(1) - 3.0935) <= tol)  
% assert(abs(probit_params.xi(52) - 0.0889) <= tol) 
% assert(abs(probit_params.indiv_lik_probit_class(1,50) - 7.3023e-04) <= tol) 
% assert(abs(probit_params.indiv_lik_probit_class(200,1) - 1) <= tol)
assert(abs(MCMC_out.pi(1,7) - 0.2076) <= tol)
assert(MCMC_out.c_i(1,1) == 22)
assert(abs(MCMC_out.theta(1,1,50,2) - 0.2108) <= tol) 
assert(abs(MCMC_out.theta(1,1,1,1) - 0.6526) <= tol)
assert(abs(MCMC_out.xi(1,1) - 3.0935) <= tol)  
assert(abs(MCMC_out.loglik(1) + 799.3670) <= tol)
assert(abs(OFMM_params.pi(35) - 0.0819) <= tol)
assert(abs(OFMM_params.theta(1,1,1) - 0.6057) <= tol)
assert(abs(OFMM_params.theta(4,50,2) - 0.1446) <= tol) 
assert(isequal(OFMM_params.c_i(1:5), [22; 35; 7; 7; 29]))
assert(OFMM_params.n_ci(7) == 38)
assert(abs(probit_params.xi(1) - 1.7924) <= tol)  
assert(abs(probit_params.xi(52) - 2.3128) <= tol) 
assert(abs(probit_params.indiv_lik_probit_class(1,50) - 2.0193e-05) <= tol) 
assert(abs(probit_params.indiv_lik_probit_class(200,1) - 1) <= tol)
assert(abs(MCMC_out.pi(9,35) - 0.0819) <= tol)
assert(MCMC_out.c_i(9,200) == 7)
assert(abs(MCMC_out.theta(9,1,50,2) - 0.5777) <= tol)
assert(abs(MCMC_out.xi(9,52) - 2.3128) <= tol)
assert(abs(MCMC_out.loglik(9) + 537.2099) <= tol)

%% Test 5: wtd_run_MCMC function
n_runs = 20;
burn = 0;
thin = 1;
MCMC_out.pi = zeros(n_runs / thin, k_max);  
MCMC_out.theta = zeros(n_runs / thin, data_vars.p, k_max, data_vars.d_max); 
MCMC_out.c_i = zeros(n_runs / thin, data_vars.n);  
MCMC_out.xi = zeros(n_runs / thin, p_cov);  
MCMC_out.z_i = zeros(n_runs / thin, data_vars.n);
MCMC_out.loglik = zeros(n_runs / thin, 1); 
[MCMC_out, OFMM_params, probit_params] = wtd_run_MCMC(data_vars, ...
        OFMM_params, probit_params, n_runs, burn, thin, k_max, q_dem, ...
        p_cov, alpha, eta, mu_0, Sig_0); 
assert(abs(MCMC_out.pi(1,7) - 0.2782) <= tol)
assert(MCMC_out.c_i(1,1) == 22)
assert(abs(MCMC_out.theta(1,1,50,2) - 0.2341) <= tol) 
assert(abs(MCMC_out.theta(1,1,1,1) - 0.8432) <= tol)
assert(abs(MCMC_out.xi(1,1) - 2.3219) <= tol)  
assert(abs(MCMC_out.loglik(1) + 505.1303) <= tol)
assert(abs(OFMM_params.pi(35) - 2.4707e-17) <= tol)
assert(abs(OFMM_params.theta(1,1,1) - 1.5570e-04) <= tol)
assert(abs(OFMM_params.theta(4,50,2) - 0.6022) <= tol) 
assert(isequal(OFMM_params.c_i(1:5), [34; 34; 14; 49; 34]))
assert(OFMM_params.n_ci(7) == 0)
assert(abs(probit_params.xi(1) - 0.3800) <= tol)  
assert(abs(probit_params.xi(52) - 0.3576) <= tol) 
assert(abs(probit_params.indiv_lik_probit_class(1,50) - 0.5049) <= tol) 
assert(abs(probit_params.indiv_lik_probit_class(200,1) - 0.6388) <= tol)
assert(abs(MCMC_out.pi(20,35) - 2.4707e-17) <= tol)
assert(MCMC_out.c_i(20,200) == 15)
assert(abs(MCMC_out.theta(20,1,50,2) - 0.7071) <= tol)
assert(abs(MCMC_out.xi(20,52) - 0.3576) <= tol)
assert(abs(MCMC_out.loglik(20) + 526.9837) <= tol)

k_fixed = round(median(sum(MCMC_out.pi > 0.05, 2)));
assert(k_fixed == 4) 

%% Test 6: Fixed sampler MCMC output
% Initialize OFMM model using fixed number of classes
clear OFMM_params probit_params MCMC_out;           % Reduce memory burden
sp_k = k_fixed;                   % Denom constant to restrict alpha size and num classes for OFMM model
alpha = ones(1, k_fixed) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
OFMM_params = init_OFMM_params(data_vars, k_fixed, alpha, eta);
% Initialize probit model using fixed number of classes
p_cov = k_fixed + S;                        % Number of covariates in probit model
mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVInvGamma(shape=5/2, scale=5/2)
Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
probit_params = init_probit_params(data_vars, k_fixed, q_dem, mu_0, Sig_0);
% Run MCMC algorithm using fixed number of classes
[MCMC_out, OFMM_params, probit_params] = wtd_run_MCMC(data_vars, ...
    OFMM_params, probit_params, n_runs, burn, thin, k_fixed, q_dem, ...
    p_cov, alpha, eta, mu_0, Sig_0);
assert(abs(MCMC_out.pi(1,3) - 1.9785e-14) <= tol)
assert(MCMC_out.c_i(1,1) == 2)
assert(abs(MCMC_out.theta(1,1,3,2) - 0.2398) <= tol) 
assert(abs(MCMC_out.theta(1,1,1,1) - 0.3878) <= tol)
assert(abs(MCMC_out.xi(1,1) - 2.0385) <= tol)  
assert(abs(MCMC_out.loglik(1) + 959.4198) <= tol)
assert(abs(OFMM_params.pi(3) - 0.4146) <= tol)
assert(abs(OFMM_params.theta(1,1,1) - 0.8030) <= tol)
assert(abs(OFMM_params.theta(4,3,2) - 0.0780) <= tol) 
assert(isequal(OFMM_params.c_i(1:5), [1; 3; 2; 2; 1]))
assert(abs(OFMM_params.n_ci(3) - 81.25) <= tol)
assert(abs(probit_params.xi(1) - 1.3748) <= tol)  
assert(abs(probit_params.xi(3) + 1.1391) <= tol) 
assert(abs(probit_params.indiv_lik_probit_class(1,3) - 0.4068) <= tol) 
assert(abs(probit_params.indiv_lik_probit_class(200,1) - 0.6400) <= tol)
assert(abs(MCMC_out.pi(20,3) - 0.4146) <= tol)
assert(MCMC_out.c_i(20,200) == 2)
assert(abs(MCMC_out.theta(20,1,3,2) - 0.0385) <= tol)
assert(abs(MCMC_out.xi(20,3) + 1.1391) <= tol)
assert(abs(MCMC_out.loglik(20) + 522.8844) <= tol)

%% Test 7: post_process function
post_MCMC_out = post_process(MCMC_out, data_vars, S);
assert(post_MCMC_out.k_med == 3)
assert(abs(post_MCMC_out.tree(1,1,1) - 137) <= tol)
assert(abs(post_MCMC_out.tree(199,1,1) - 394) <= tol)
assert(abs(post_MCMC_out.pi(1,3) -  0.3333) <= tol)  
assert(abs(post_MCMC_out.pi(20,3) -  0.2359) <= tol) 
assert(abs(post_MCMC_out.theta(20,1,1,1) - 0.9615) <= tol)
assert(abs(post_MCMC_out.theta(1,3,1,2) - 0.3152) <= tol)
assert(abs(post_MCMC_out.xi(1,1) - 2.0385) <= tol)
assert(abs(post_MCMC_out.xi(20,3) + 0.9999) <= tol) 
assert(abs(post_MCMC_out.loglik(20) + 522.8844) <= tol)

%% Test 8: analyze_results function
analysis = analyze_results(MCMC_out, post_MCMC_out, data_vars, q_dem, S, p_cov);
assert(abs(analysis.pi_med(1) - 0.5910) <= tol)
assert(abs(analysis.theta_med(1,1,1) - 0.8945) <= tol)
assert(abs(analysis.theta_med(4,2,2) - 0.7828) <= tol)
assert(abs(analysis.xi_med(1) - 1.1490) <= tol)
assert(abs(analysis.xi_med(4) + 0.1492) <= tol)
assert(analysis.k_red == 2)
assert(isequal(analysis.c_i(1:5), [1; 1; 2; 2; 1]))
assert(abs(analysis.max_prob_ci(1) - 0.9986) <= tol)
assert(abs(analysis.max_prob_ci(200) - 0.9764) <= tol)
assert(abs(analysis.Phi_med(1) - 0.7641) <= tol)
assert(abs(analysis.Phi_med(200) - 0.8413) <= tol)
assert(abs(sum(analysis.loglik_med) + 498.8875) <= tol)  
analysis.dic6
analysis.aebic


