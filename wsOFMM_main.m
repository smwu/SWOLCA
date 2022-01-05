
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted Supervised OFMM Main
% Programmer: SW    
% Data: Simulations       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% Load simulated data
sim_n = 1;   % Simulation number
scenario = 4;  % Weighting scenario
samp_data = importdata(strcat('simdata_wsRPC_scen', num2str(scenario), '_iter',num2str(sim_n),'.mat'));

%% Get data variables
k_max = 50;    % Upper limit for number of classes
data_vars = get_data_vars(samp_data);

%% Initialize priors and variables for OFMM model
sp_k = 50;                      % Denom constant to restrict alpha size and num classes for OFMM model
alpha = ones(1, k_max) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
eta = ones(1, data_vars.d_max); % Hyperparam for item response probs. 
OFMM_params = init_OFMM_params(data_vars, k_max, alpha, eta);

%% Initialize priors and variables for probit model
q_dem = [];                                 % Matrix of demographic covariates in cell-means format. Default is empty.
S = size(q_dem, 2);                         % Number of demographic covariates in the probit model
p_cov = k_max + S;                          % Number of covariates in probit model
mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVGamma(shape=5/2, scale=5/2)
Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 
probit_params = init_probit_params(data_vars, k_max, q_dem, mu_0, Sig_0);

%% Run MCMC algorithm to obtain posteriors and save output
% n_runs = 25000;  % Number of MCMC iterations
% burn = 15000;    % Burn-in period
% thin = 5;        % Thinning factor
n_runs = 50;  % Number of MCMC iterations
burn = 30;    % Burn-in period
thin = 2;        % Thinning factor
[MCMC_out, OFMM_params, probit_params] = run_MCMC(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_max, q_dem, p_cov, alpha, eta, mu_0, Sig_0, S, sim_n, scenario);

%% Post-processing to recalibrate labels and remove extraneous empty classes
post_MCMC_out = post_process(MCMC_out, data_vars, S);

%% Obtain posterior estimates, reduce number of classes, analyze results, and save output
analysis = analyze_results(post_MCMC_out, data_vars, q_dem, S, sim_n, scenario);


%% LOCAL FUNCTIONS

% get_data_vars takes in sample data and outputs relevant variables
% Inputs: samp_data structural array with at least the following columns:
%   sample_data: food intake data as a matrix
%   true_y: outcome data as a vector
%   sample_wt: survey weights as a vector
% Outputs: data_vars structural array with the following fields:
%   food: matrix of food intake data
%   n: number of individuals
%   p: number of food items
%   d_max: max number of consumption levels over all items
%   d: vector of max number of levels for each food item
%   y: vector of outcomes
%   wt_kappa: vector of normalized weights
%   wt_kappa_mat: matrix of normalized weights, replicated across items
%   lin_idx: vector of linear indices for unique item-response combos
function data_vars = get_data_vars(samp_data)
    data_vars.food = samp_data.sample_data;
    [data_vars.n, data_vars.p] = size(data_vars.food);
    data_vars.d_max = max(data_vars.food(:));    % Max number of levels over all items
    data_vars.d = max(data_vars.food);           % Max number of levels for each food item. 
    data_vars.y = samp_data.true_y;            % Outcome
    
    % Normalized weights
    kappa = sum(samp_data.sample_wt) / data_vars.n;    % Normalization constant. If weights sum to N, this is 1/(sampl frac)
    data_vars.wt_kappa = samp_data.sample_wt / kappa;  % Normalized weight to sum to n instead of N
    data_vars.wt_kappa_mat = repmat(data_vars.wt_kappa, [1, data_vars.p]); % nxp weight matrix. Weights replicated across food items
    
    % Item-response combinations
    idz = repmat(1:data_vars.p, data_vars.n, 1);  % nxp matrix of item id's. Replicates 1:p row n times
    idz = idz(:);                                 % cols of idz, concat by col into vector of length np. 1(x n),2(x n),...
    x_d = data_vars.food(:);                      % Food item responses as a long vector (corresponds to idz) 
    % lin_idx: np vector of unique item-response (idz x y_d) combo indices
    % Each unique combo has corresp linear index. Those with same combo have same lin_idx value
    data_vars.lin_idx = sub2ind([data_vars.p, data_vars.d_max], idz, x_d);
end


% init_OFMM_params initializes priors and variables for the OFMM model given hyperparameters
% Inputs: 
%   data_vars: output from get_data_vars function
%   k_max: max number of latent classes
%   alpha: vector of hyperparams for class membership probs, pi. Default is alpha=1/50
%   eta: vector of hyperparams for item response probs, theta. Default is eta=1
% Outputs: OFMM_params structural array with the following fields:
%   pi: vector prior for class membership probs. Dist: Dir(alpha)
%   theta: vector prior for item response probabilities. Dist: Dir(eta) per class and item
%   c_i: vector of initial class assignments (randomly generated)
%   n_ci: vector of num indivs initially assigned to each cluster
function OFMM_params = init_OFMM_params(data_vars, k_max, alpha, eta)
    % Prior for class membership probs
    OFMM_params.pi = drchrnd(alpha, 1); 
    
    % Prior for item response probabilities
    OFMM_params.theta = zeros(data_vars.p, k_max, data_vars.d_max);  % Initialize array of item response probs
    for k = 1:k_max                % For each class
        for j = 1:data_vars.p      % For each food item
            d_j = data_vars.d(j);  % Max number of levels for food item j
            OFMM_params.theta(j, k, 1:d_j) = drchrnd(eta(1:d_j), 1); % For each class and item, draw from Dir(eta)
        end
    end
    
    % Random initialization of class assignments
    rr = unifrnd(0, 1, [data_vars.n, 1]);     % Vector of Unif(0,1) draws
    pisums = [0 cumsum(OFMM_params.pi)];      % Cum sum of class membership probs, from 0 to 1
    OFMM_params.c_i = zeros(data_vars.n, 1);  % Initialize vector of class assignments
    x_ci = zeros(data_vars.n, k_max);         % Initialize matrix of cell-means class assignments 
    for k = 1:k_max
        ind = rr > pisums(k) & rr <= pisums(k + 1);  % T/F indicating those w/ rr in the kth pisums interval
        OFMM_params.c_i(ind == 1) = k;               % Class assignment for those w/ ind == TRUE is set to k
        x_ci(:, k) = ind;                            % For each indiv, 1 in col corresp to class assignment
    end 
    OFMM_params.n_ci=sum(x_ci);              %Vector of num indivs assigned to each cluster
end

% init_probit_params initializes priors and variables for the probit model given hyperparameters
% Inputs: 
%   data_vars: output from get_data_vars function
%   k_max: max number of latent classes
%   q_dem: matrix of demographic covariates to include in regression, in cell-means format
%   mu_0: vector of mean hyperparams for probit coefficients, xi. 
%         Default is normrnd(0,1,[p_cov,1]). MVN(0,1) with indep components
%   Sig_0: matrix of variance hyperparams for probit coefficients, xi.
%          Default is diag(1./gamrnd(5/2,2/5,[pcov,1])). Gamma(shape=5/2,scale=5/2) with indep components
% Outputs: probit_params structural array with the following fields:
%   xi: vector prior for probit coefficients. Dist: MVN(mu_0, sig_0)
%   indiv_lik_probit: matrix of indiv probit lik components used for posterior class membership
function probit_params = init_probit_params(data_vars, k_max, q_dem, mu_0, Sig_0)
    % Prior for model params, drawn from MVN(mu0, Sig0)
    probit_params.xi = mvnrnd(mu_0, Sig_0);  
    
    % Probit likelihood component of posterior class membership P(c_i|-)
    probit_params.indiv_lik_probit = zeros(data_vars.n, k_max);    % Initialize matrix for individual probit likelihood contribution
    q_class = zeros(data_vars.n, k_max);                           % Initialize design matrix for class membership
    for k = 1:k_max
       q_class(:, k) = ones(data_vars.n, 1);                       % Temporarily assign all indivs to class k
       Q_temp = [q_dem q_class];                                   % Form temporary covariate design matrix
       Phi_temp = normcdf(Q_temp * transpose(probit_params.xi));   % Vector of P(y_i=1|Q)
       Phi_temp(Phi_temp == 1) = 1 - 1e-10;                        % Adjust extremes
       Phi_temp(Phi_temp == 0) = 1e-10;
       loglik_probit = data_vars.y .* log(Phi_temp) + (1-data_vars.y) .* log(1-Phi_temp); % Log-lik helps avoid rounding issues
       probit_params.indiv_lik_probit(:, k) = exp(loglik_probit);  % Probit likelihood for assigning indivs to class k
    end
end

% MCMC_out runs the Gibbs Sampler MCMC algorithm to obtain posteriors
% Inputs: 
%   data_vars: output from get_data_vars function
%   OFMM_params: output from init_OFMM_params function
%   probit_params: output from init_probit_params function
%   n_runs: number of MCMC iterations
%   burn: burn-in period
%   thin: thinning factor
%   k_max: max number of latent classes
%   q_dem: matrix of demographic covariates to include in regression, in cell-means format
%   p_cov: number of covariates in probit model
%   alpha: vector hyperparam for class membership probs, pi
%   eta: vector hyperparam for item-response probs, theta
%   mu_0: vector mean hyperparam for probit coefficients, xi 
%   Sig_0: matrix variance hyperparam for probit coefficients, xi
%   S: number of demographic covariates in the probit model
%   sim_n: simulation number indicator
%   scenario: weighting scenario indicator
% Outputs: returns and saves MCMC_out structural array with the following fields:
%   pi: matrix of class membership probabilities 
%   theta: 4-D array of item-response probabilities
%   c_i: matrix of class assignments
%   xi: matrix of probit model coefficients
%   z_i: matrix of latent probit outcomes
%   loglik: vector of log likelihoods
function [MCMC_out, OFMM_params, probit_params] = run_MCMC(data_vars, OFMM_params, probit_params, n_runs, burn, thin, k_max, q_dem, p_cov, alpha, eta, mu_0, Sig_0, S, sim_n, scenario)
    % Initialize output variables
    MCMC_out.pi = zeros(n_runs / thin, k_max);  
    MCMC_out.theta = zeros(n_runs / thin, data_vars.p, k_max, data_vars.d_max); 
    MCMC_out.c_i = zeros(n_runs / thin, data_vars.n);  
    MCMC_out.xi = zeros(n_runs / thin, p_cov);  
    MCMC_out.z_i = zeros(n_runs / thin, data_vars.n);
    MCMC_out.loglik = zeros(n_runs / thin, 1); 
    
    % Run MCMC iterations
    tic                  % Start timer
    for iter = 1:n_runs  % For each MCMC run
        % Update parameters and variables
        [MCMC_out, OFMM_params, probit_params] = update_MCMC(MCMC_out, data_vars, OFMM_params, probit_params, thin, k_max, q_dem, alpha, eta, mu_0, Sig_0, iter);
    
        % Relabel classes every 10 iterations to encourage mixing
        if mod(iter, 10) == 0
            new_order = randperm(k_max);                               % New labels for latent classes
            new_ci = OFMM_params.c_i;
            for k = 1:k_max
                new_ci(OFMM_params.c_i == k) = new_order(k);           % Class assignments with new labels
            end
            OFMM_params.c_i = new_ci;                                  % Relabel class assignments
            OFMM_params.theta = OFMM_params.theta(:, new_order, :);    % Relabel item-response probabilities
            probit_params.indiv_lik_probit = probit_params.indiv_lik_probit(:, new_order);  % Relabel indiv probit likelihood
            OFMM_params.pi = OFMM_params.pi(new_order);                % Relabel class probabilities
            probit_params.xi((S+1):end) = probit_params.xi(S + new_order); % Relabel probit model coefficients
            alpha = alpha(new_order);                                  % Relabel hyperparam for pi
            mu_0((S+1):end) = mu_0(S + new_order);                     % Relabel mean hyperparam for xi
            idx = find(Sig_0);                                         % Diagonal elements of Sig_0
            Sig_0(idx((S+1):end)) = Sig_0(idx(S + new_order));         % Relabel var hyperparam for xi
        end
    end
    MCMC_out.runtime = toc;  % Stop timer
    
    % Discard burn-in
    MCMC_out.pi = MCMC_out.pi((burn / thin) + 1:end, :);
    MCMC_out.theta = MCMC_out.theta((burn / thin) + 1:end, :, :, :);
    MCMC_out.c_i = MCMC_out.c_i((burn / thin) + 1:end, :);
    MCMC_out.xi = MCMC_out.xi((burn / thin) + 1:end, :);  
    MCMC_out.z_i = MCMC_out.z_i((burn / thin) + 1:end);
    MCMC_out.loglik = MCMC_out.loglik((burn / thin) + 1:end); 
    
    % Save output
    save(strcat('wsOFMM_MCMC_out','_',num2str(scenario),'_',num2str(sim_n)), 'MCMC_out');
end

% MCMC_update updates the posterior distributions of the parameters and variables
% Inputs: 
%   MCMC_out: output initializes in MCMC_run
%   data_vars: output from get_data_vars function
%   OFMM_params: output from init_OFMM_params function
%   probit_params: output from init_probit_params function
%   thin: thinning factor
%   k_max: max number of latent classes
%   q_dem: matrix of demographic covariates to include in regression, in cell-means format
%   alpha: vector hyperparam for class membership probs, pi
%   eta: vector hyperparam for item-response probs, theta
%   mu_0: vector mean hyperparam for probit coefficients, xi 
%   Sig_0: matrix variance hyperparam for probit coefficients, xi
%   iter: MCMC iteration indicator
% Outputs: returns updated OFMM_params and probit_params structural arrays,
% as well as updated MCMC_out structural array with the following fields:
%   pi: matrix of class membership probabilities 
%   theta: 4-D array of item-response probabilities
%   c_i: matrix of class assignments
%   xi: matrix of probit model coefficients
%   z_i: matrix of latent probit outcomes
%   loglik: vector of log likelihoods
function [MCMC_out, OFMM_params, probit_params] = update_MCMC(MCMC_out, data_vars, OFMM_params, probit_params, thin, k_max, q_dem, alpha, eta, mu_0, Sig_0, iter) 
    % Update posterior class membership probs, pi, incorporating normalized weights
    for k = 1:k_max  % For each latent class
        % Get updated weighted num indivs assigned to class k
        OFMM_params.n_ci(k) = sum(data_vars.wt_kappa(OFMM_params.c_i == k));  
    end
    alpha_post = alpha + OFMM_params.n_ci;    % Posterior hyperparam is alpha(k) + n_ci(k)
    OFMM_params.pi = drchrnd(alpha_post, 1);  % Draw from posterior dist: Dir(alpha_post)
    
    % Update class membership, c_i
 	prob_ci_numer = zeros(data_vars.n, k_max);  % Initialize matrix of P(c_i=k|-) numerator for each indiv and class
    for k = 1:k_max
        % Matrix of unique item-level probs for class k
        theta_k = reshape(OFMM_params.theta(:, k, :), data_vars.p, data_vars.d_max); 
        % Matrix of categ likelihood \prod_{r=1}^d \theta_{jr|k}^{I(x_ij=r)} for each indiv and item 
        categ_lik = reshape(theta_k(data_vars.lin_idx), [data_vars.n, data_vars.p]); 
        % Update numerator for posterior P(c_i=k|-) 
        prob_ci_numer(:, k) = OFMM_params.pi(k) * prod(categ_lik, 2) .* probit_params.indiv_lik_probit(:, k);  
    end
    % Matrix of updated posterior P(c_i=k|-) for each indiv and class
    prob_ci = bsxfun(@times, prob_ci_numer, 1 ./ (sum(prob_ci_numer, 2)));
    indiv_loglik = log(sum(prob_ci_numer, 2));  % Individual log likelihood
    x_ci = mnrnd(1, prob_ci);                   % Matrix of updated c_i drawn from Mult(1, prob_ci). Each row is draw for an indiv
    [row, col] = find(x_ci);                    % Row and col indices of each nonzero element of x_ci. col is class assignment for each indiv
    x_gc = [row col]; x_gc = sortrows(x_gc, 1); % Matrix with 1st and 2nd col the row and col indices
    OFMM_params.c_i = x_gc(:, 2);               % Vector of updated class assignment, c_i, for each individual
    
    % Update posterior item response probabilities, theta, incorporating normalized weights    
    n_theta = zeros(data_vars.p, data_vars.d_max); % Initialize matrix of num indivs with each item response level
    for k = 1:k_max  
        c_i_mat = repmat(OFMM_params.c_i, [1, data_vars.p]);  % Rep c_i's p times to form matrix of classes for each indiv and item
        class_wt = (c_i_mat == k) .* data_vars.wt_kappa_mat;  % Matrix of normalized weights for indivs assigned to class k. All other rows set to 0
        for r = 1:data_vars.d_max                             % For each consumption level r
            % Vector of num indivs with level r and class k, for each item (i.e., sum I(x_ij=r, c_i=k) over indivs)
            n_theta(:, r) = sum((data_vars.food == r) .* class_wt)'; 
        end
        for j = 1:data_vars.p                                               % For each food item
            dj = data_vars.d(j);                                  % Num response levels for item j
            eta_post = eta(1:dj) + n_theta(j, 1:dj);              % Posterior hyperparam is eta[r] + n_theta(j,r)
            OFMM_params.theta(j, k, 1:dj) = drchrnd(eta_post, 1); % Draw from posterior dist: Dir(eta), for every item and class
        end
    end        
    
    % Update latent probit variable, z_i
    z_i = zeros(data_vars.n, 1);                 % Initialize latent probit variable, z_i
    Q = [q_dem x_ci];                            % Design matrix with demographic covariates and updated latent classes
    lin_pred = Q * transpose(probit_params.xi);  % Linear predictor, Q*xi. Mean of truncated normal dist
    % For cases, z_i ~ TruncNormal(mean=Q*xi, var=1, low=0, high=Inf)
    z_i(data_vars.y == 1) = truncnormrnd(1, lin_pred(data_vars.y == 1), 1, 0, Inf); 
    % For controls, z_i ~ TruncNormal(mean=Q*xi, var=1, low=-Inf, high=0)
    z_i(data_vars.y == 0) = truncnormrnd(1, lin_pred(data_vars.y == 0), 1, -Inf, 0); 
    if sum(z_i == Inf) > 0                       % Control extremes
        z_i(z_i == Inf) = norminv(1 - 1e-10);
    end 
    if sum(z_i == -Inf) > 0
        z_i(z_i == -Inf) = norminv(1e-10);
    end    
    
    % Update posterior probit model coefficients, xi, incorporating normalized weights    
    W_tilde = diag(data_vars.wt_kappa);                              % Diagonal normalized weight matrix
    Sig_post = inv(Sig_0) + (transpose(Q)*W_tilde*Q);                % Precision of posterior Normal dist for xi
    mu_right = (Sig_0*mu_0) + (transpose(Q)*W_tilde*z_i); % Right part of mean of posterior Normal dist for xi
    mu_post = Sig_post \ mu_right;                                   % Mean of posterior Normal dist for xi
    probit_params.xi = mvnrnd(mu_post, inv(Sig_post));               % Draw from posterior dist: N(mu_post, Sig_post^{-1})    
    
    % Update probit likelihood component of posterior class membership P(c_i|-)
    q_class = zeros(data_vars.n, k_max);                           % Initialize design matrix for class membership
    for k = 1:k_max
       q_class(:, k) = ones(data_vars.n, 1);                       % Temporarily assign all indivs to class k
       Q_temp = [q_dem q_class];                                   % Form temporary covariate design matrix
       Phi_temp = normcdf(Q_temp * transpose(probit_params.xi));   % Vector of P(y_i=1|Q)
       Phi_temp(Phi_temp == 1) = 1 - 1e-10;                        % Adjust extremes
       Phi_temp(Phi_temp == 0) = 1e-10;
       % Update indiv probit log-lik for class k. Helps avoid rounding issues
       loglik_probit = data_vars.y .* log(Phi_temp) + (1-data_vars.y) .* log(1-Phi_temp);
       probit_params.indiv_lik_probit(:, k) = exp(loglik_probit);  % Update indiv probit likelihood for class k
    end    
    
    % Store posterior values if the iteration is selected based on thinning  
    if mod(iter, thin) == 0
        MCMC_out.pi(iter / thin, :) = OFMM_params.pi;
        MCMC_out.c_i(iter / thin, :) = OFMM_params.c_i;
        MCMC_out.theta(iter / thin, :, :, :) = OFMM_params.theta; 
        MCMC_out.loglik(iter / thin) = sum(indiv_loglik); 
        MCMC_out.z_i(iter / thin, :) = z_i;
        MCMC_out.xi(iter / thin, :) = probit_params.xi; 
    end
end


% post_process conducts post-processing to remove extraneous empty classes and relabel class assignments
% Inputs: 
%   MCMC_out: output initialized in MCMC_run
%   data_vars: output from get_data_vars function
%   S: number of demographic covariates in the probit model 
% Outputs: returns post_MCMC_out structural array with the following fields:
%   k_med: Median num classes w/ more than >5% indivs, calculated over MCMC iterations
%   tree: Hierarchical clustering dendrogram tree using distance criterion to find num classes
%   pi: Matrix of reordered class membership probs
%   theta: Matrix reordered item-response probs
%   xi: Matrix of reordered probit model coefficients
function post_MCMC_out = post_process(MCMC_out, data_vars, S) 
    % Remove extraneous empty classes and relabel class assignments  
    m = size(MCMC_out.pi, 1);                                 % Num stored output iterations
    post_MCMC_out.k_med = median(sum(MCMC_out.pi > 0.05, 2)); % Median num classes w/ more than >5% indivs, over all iterations
    pd = pdist(transpose(MCMC_out.c_i), 'hamming');           % Pairwise distance of num times two indivs were in same class together
    cdiff = squareform(pd);                                   % Convert pd into square matrix
    post_MCMC_out.tree = linkage(cdiff,'complete');           % Hierarchical clustering dendrogram tree using distance criterion to find num classes
    red_ci = cluster(post_MCMC_out.tree, 'maxclust', post_MCMC_out.k_med); % Vector of new class assignments for each indiv
    relabel_ci = zeros(m, post_MCMC_out.k_med);                            % Initialize matrix of relabeled class assignments after removing empty classes
    for k = 1:post_MCMC_out.k_med                                          % For each non-empty class
        % For each MC run, among indivs assigned to new class k, which original class was most common 
        relabel_ci(:, k) = mode(MCMC_out.c_i(:, red_ci == k), 2); 
    end
    
    % Reorder parameter estimates to match the correct new class assignment
    post_MCMC_out.pi = zeros(m, post_MCMC_out.k_med);         % Initialize parameter estimates corresponding to the non-empty classes
    post_MCMC_out.theta = zeros(m, data_vars.p, post_MCMC_out.k_med, data_vars.d_max);
    post_MCMC_out.xi = zeros(m, S + post_MCMC_out.k_med);
    for iter = 1:m                                            % For each stored MC run
        iter_order = relabel_ci(iter, :);                     % Vector of new class assignments (mode orig classes after removing empty classes)
        pi_order = MCMC_out.pi(iter, iter_order);             % Update class probs corresp to new class assignments
        post_MCMC_out.pi(iter, :) = pi_order / sum(pi_order); % Renormalize to sum to 1
        % Update item-response probs corresp to new class assignments
        post_MCMC_out.theta(iter, :, :, :) = MCMC_out.theta(iter, :, iter_order, :);  
        % Update probit model coefs corresp to new class assignments. First few coefs corresp to dem vars and are not changed
        post_MCMC_out.xi(iter, :) = [MCMC_out.xi(iter, 1:S) MCMC_out.xi(iter, S + iter_order)];  
    end
end


% analyze_results obtains reduced number of classes, posterior estimates of the params, and analysis metrics
% Inputs: 
%   post_MCMC_out: output from post_process function
%   data_vars: output from get_data_vars function
%   q_dem: matrix of demographic covariates to include in regression, in cell-means format
%   S: number of demographic covariates in the probit model
%   sim_n: simulation number indicator
%   scenario: weighting scenario indicator
% Outputs: returns analysis structural array with the following fields:
%   pi_med: Vector of posterior median estimate across MCMC runs for class membership probs
%   theta_med: Array of posterior median estimate across MCMC runs for item-response probs
%   xi_med: Vector of posterior median estimate across MCMC runs for probit model coefs
%   k_red: Reduced number of classes after removing duplicates
%   indiv_lik_probit: Matrix of indiv probit likelihood using unique classes and estimate for xi
%   loglik: Log likelihood
%   c_i: Vector of class assignment for each individual using posterior median estimates
%   Phi_med: estimated probit model mean using posterior median estiates 
%   loglik_mean: Log-lik using posterior median estimates of all params
%   dic: Assesses model goodness-of-fit
function analysis = analyze_results(post_MCMC_out, data_vars, q_dem, S, sim_n, scenario)
    % Calculate posterior median estimates across MCMC runs
    analysis.pi_med = median(post_MCMC_out.pi);  % Vector of median class probs 
    % Array of median item-response probs
    analysis.theta_med = reshape(median(post_MCMC_out.theta), [data_vars.p, post_MCMC_out.k_med, data_vars.d_max]);  
    analysis.xi_med = median(post_MCMC_out.xi);  % Vector of median probit coefs

    % Identify unique classes    
    % Matrix of most likely consump level based on item-response probs, for each item and class across MCMC runs
    [~, ind0] = max(analysis.theta_med, [], 3); 
    t_ind0 = transpose(ind0);
    [~, ia, ~] = unique(t_ind0, 'rows');  % Vector of class indices corresponding to unique classes
    analysis.k_red = length(ia);          % Number of unique classes
    
    % Posterior median estimates for the unique classes
    analysis.pi_med = analysis.pi_med(ia) / sum(analysis.pi_med(ia)); % Class probs for the unique classes, normalized 
    % Item-response probs for the unique classes, normalized to sum to 1 across consumption levels
    analysis.theta_med = analysis.theta_med(:, ia, :); 
    analysis.theta_med = analysis.theta_med ./ sum(analysis.theta_med, 3);  
    analysis.xi_med = [analysis.xi_med(1:S) analysis.xi_med(S + ia)]; % Probit coefs for the unique classes
    
    % Update probit likelihood using unique classes and updated posterior median estimates for xi
    analysis.indiv_lik_probit = zeros(data_vars.n, analysis.k_red); % Initialize indiv probit lik matrix for each indiv and class
    q_class = zeros(data_vars.n, analysis.k_red);               % Initialize design matrix for class membership
    for k = 1:analysis.k_red
       q_class(:, k) = ones(data_vars.n, 1);                    % Temporarily assign all indivs to class k
       Q_temp = [q_dem q_class];                                % Form temporary covariate design matrix
       Phi_temp = normcdf(Q_temp * transpose(analysis.xi_med)); % Vector of P(y_i=1|Q)
       Phi_temp(Phi_temp == 1) = 1 - 1e-10;                     % Adjust extremes
       Phi_temp(Phi_temp == 0) = 1e-10;
       % Update indiv probit log-lik for class k. Helps avoid rounding issues
       loglik_probit = data_vars.y .* log(Phi_temp) + (1-data_vars.y) .* log(1-Phi_temp);
       analysis.indiv_lik_probit(:, k) = exp(loglik_probit);    % Update indiv probit likelihood for class k
    end  
    
    % Update class membership, c_i, using unique classes and updated posterior median estimates for pi and theta
 	prob_ci_numer = zeros(data_vars.n, analysis.k_red);  % Initialize indiv OFMM lik matrix for each indiv and class
    for k = 1:analysis.k_red
        % Matrix of unique item-level probs for class k
        theta_k = reshape(analysis.theta_med(:, k, :), data_vars.p, data_vars.d_max); 
        % Matrix of categ likelihood \prod_{r=1}^d \theta_{jr|k}^{I(x_ij=r)} for each indiv and item 
        categ_lik = reshape(theta_k(data_vars.lin_idx), [data_vars.n, data_vars.p]); 
        % Update numerator for posterior P(c_i=k|-) (indiv likelihood)
        prob_ci_numer(:, k) = analysis.pi_med(k) * prod(categ_lik, 2) .* analysis.indiv_lik_probit(:, k);   
    end
    % Matrix of updated posterior P(c_i=k|-) for each indiv and class
    prob_ci = bsxfun(@times, prob_ci_numer, 1 ./ (sum(prob_ci_numer, 2)));
    analysis.loglik = sum(log(sum(prob_ci_numer, 2)));  % Log likelihood
    x_ci = mnrnd(1, prob_ci);                   % Matrix of updated c_i drawn from Mult(1, prob_ci). Each row is draw for an indiv
    [row, col] = find(x_ci);                    % Row and col indices of each nonzero element of x_ci. col is class assignment for each indiv
    x_gc = [row col]; x_gc = sortrows(x_gc, 1); % Matrix with 1st and 2nd col the row and col indices
    analysis.c_i = x_gc(:, 2);                  % Vector of updated class assignment, c_i, for each individual
    
    % Obtain probit model mean using median parameter estimates
    Q_med = [q_dem x_ci];                             % Design matrix with updated class memberships using median estimates   
    analysis.Phi_med = normcdf(Q_med * transpose(analysis.xi_med)); % Probit model mean using median estimates
    analysis.Phi_med(analysis.Phi_med == 0) = 1e-10;  % Adjust extremes
    analysis.Phi_med(analysis.Phi_med == 1) = 1 - 1e-10;
        
    % Calculate log-lik using posterior median estimates of all params
    loglik_probit_mean = data_vars.y .* log(analysis.Phi_med) + (1 - data_vars.y) .* log(1 - analysis.Phi_med);
    analysis.loglik_mean = analysis.pi_med(k) * prod(categ_lik, 2) .* loglik_probit_mean;
    
    % DIC is a metric for model goodness of fit. It compares two values:
        % 1) median(loglik): calculate log-lik at each MCMC iter then get the median
        % 2) loglik_mean: get posterior estimates of each param, then use these to calculate the log-lik
    analysis.dic = -4 * median(analysis.loglik) + 2 * sum(analysis.loglik_mean);    
    
    % Save parameter estimates and analysis results
    save(strcat('wsOFMM_Results','_',num2str(scenario),'_',num2str(sim_n)), 'post_MCMC_out', 'analysis');
end

% drchrnd generates a random sample from a Dirichlet distribution
% Inputs:
%   a: Vector of hyperparameter values
%   n: Size of the sample to draw
% Outputs:
%   r: Random sample vector of size n drawn from a Dir(a) distribution
function r = drchrnd(a,n)
    p = length(a);
    r = gamrnd(repmat(a, n, 1), 1, n, p);
    r = r ./ repmat(sum(r, 2), 1, p);
end

% truncnormrnd generates a random sample from a truncated bormal
% distribution. The statistics toolbox is used.
% Inputs:
%   N: Size of the sample to draw
%   mu: Mean of the underlying normal distribution
%   sig: Standard deviation of the underlying normal distribution
%   xlo: Low truncation point, if any
%   xhigh: High truncation point, if any
% Outputs:
%   z: Random sample vector of size N drawn from a TN(mu, sig, xlo, xhi) distribution
function z = truncnormrnd(N, mu, sig, xlo, xhi)
    if (nargin < 2) | isempty(mu)      % Default mean
      mu = 0;
    end
    if (nargin < 3) | isempty(sig)     % Default standard deviation
      sig = 1;
    end
    if (nargin < 4) | isempty(xlo)     % Default low truncation point
      xlo = -inf;
      plo = 0;
    else
      plo = normcdf((xlo - mu) / sig); % Standard normal percentile for low truncation point
    end
    if (nargin < 5) | isempty(xhi)     % Default high truncation point
      xhi = inf;
      phi = 1;
    else
      phi = normcdf((xhi - mu) / sig); % Standard normal percentile for high truncation point
    end

    if xlo > xhi  % Test if trunation points are reversed
      error 'Must have xlo <= xhi if both provided'
    end

    r = rand(N);                % Generate Unif[0,1] random deviates
    r = plo + (phi - plo) * r;  % Scale to [plo, phi]
    z = norminv(r);             % Invert through standard normal
    z = mu + z * sig;           % Apply shift and scale
end

%     % Extra code to create and save figures
%     figure; %check mixing of pi parameter
%     plot(MCMC_out.pi)
%     saveas(gcf,'wsOFMM_pis.png')
%     
%     figure; %save dendrogram of hierarchical clustering
%     dendrogram(post_MCMC_out.tree);
%     saveas(gcf,'wsOFMM_dendrogram.png')
%
%     % Compute the MSE comparing the estimated probit model to the true model 
%     y_mse = immse(analysis.Phi_med, samp_data.true_Phi)  % Phi_true must be pulled from simulated dataset
    