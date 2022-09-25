%%%%%%%%%%%%%%%%%%%%%%%%
% Supervised OFMM Main %
% Programmer: SW       %
% Data: Simulations    %
%%%%%%%%%%%%%%%%%%%%%%%%

% Outputs: saves MCMC output in a file named 'sOFMM_MCMC...' and saves 
% posterior results in a file named 'sOFMM_results...'

%% PARAMETER SETUP
scenario = 7; %% SRS with n=400
sim_n = 1; % Iter = 1
samp_n = 1; % Sample number = 1
rng(sim_n, 'twister');  % set seed

%% Load simulated data and check if file already exists
% Input directory
in_dir = strcat(pwd, '/');
% Output directory 
out_dir = in_dir;   
if samp_n > 0   % If working with sample 
    samp_data = importdata(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '_samp', num2str(samp_n), '.mat'));
else            % If working with population
    samp_data = importdata(strcat(in_dir, 'simdata_scen', num2str(scenario), '_iter', num2str(sim_n), '.mat'));
end


%% FIXED SAMPLER WITH FACTOR VARIABLE CODING

% Get data variables
data_vars.food = samp_data.X_data;
[data_vars.n, data_vars.p] = size(data_vars.food);
data_vars.d_max = max(data_vars.food(:));  % Max number of levels over all items
data_vars.d = max(data_vars.food);         % Max number of levels for each food item. 
data_vars.y = samp_data.Y_data;            % Outcome
% Normalized weights
kappa = sum(samp_data.sample_wt) / data_vars.n;    % Normalization constant. If weights sum to N, this is 1/(sampl frac)
data_vars.wt_kappa = samp_data.sample_wt / kappa;  % Normalized weight to sum to n instead of N
data_vars.wt_kappa_mat = repmat(data_vars.wt_kappa, [1, data_vars.p]); % nxp weight matrix. Weights replicated across food items
% Item-response combinations with assumed class
item_idx = repmat(1:data_vars.p, data_vars.n, 1);  % nxp matrix of item id's. Replicates 1:p row n times
item_idx = item_idx(:);                            % cols of idz, concat by col into vector of length np. 1(x n),2(x n),...
x_idx = data_vars.food(:);                         % Food item responses as a long vector (corresponds to idz) 
% lin_idx: np vector of unique item-response (idz x y_d) combo indices
% Each unique combo has corresp linear index. Those with same combo have same lin_idx value
data_vars.lin_idx = sub2ind([data_vars.p, data_vars.d_max], item_idx, x_idx);

k_fixed = 2;
k_max = k_fixed;  % Running fixed sampler
eta = ones(1, data_vars.d_max); % Hyperparam for item response probs. 
% Factor variable version
q_dem = dummyvar(samp_data.true_Si);        % Matrix of demographic covariates in cell-means format. Default contains subpop. 
S = size(q_dem, 2);                         % Number of demographic covariates in the probit model


% Initialize OFMM model using fixed number of classes 
sp_k = k_fixed;                   % Denom constant to restrict alpha size and num classes for OFMM model
alpha = ones(1, k_fixed) / sp_k;  % Hyperparam for class membership probs. R package 'Zmix' allows diff options
% Prior for class membership probs
OFMM_params.pi = drchrnd(alpha, 1); 
% Random initialization of class assignments
OFMM_params.x_ci = mnrnd(1, OFMM_params.pi, data_vars.n); % Matrix of c_i drawn from Mult(1, pi). Each row is draw for an indiv
[row, col] = find(OFMM_params.x_ci);          % Row and col indices of nonzero elements of x_ci; col is class assign for each indiv
sorted = sortrows([row col], 1);              % Sort indices in ascending row index order, to match subject index
OFMM_params.c_i = sorted(:, 2);               % Vector of class assignment, c_i, for each individual
OFMM_params.n_ci = sum(OFMM_params.x_ci);     % Vector of num indivs assigned to each cluster       
% Prior for item response probabilities
OFMM_params.theta = zeros(data_vars.p, k_max, data_vars.d_max);  % Initialize array of item response probs
for k = 1:k_max                % For each class
    for j = 1:data_vars.p      % For each food item
        d_j = data_vars.d(j);  % Max number of levels for food item j
        OFMM_params.theta(j, k, 1:d_j) = drchrnd(eta(1:d_j), 1); % For each class and item, draw from Dir(eta)
    end
end     


% Initialize probit model using fixed number of classes
% Factor variable version
p_cov = S * k_max;                        % Interactions
mu_0 = normrnd(0, 1, [p_cov, 1]);           % Mean hyperparam drawn from MVN(0,1)
Sig_0 = 1 ./ gamrnd(5/2, 2/5, [p_cov, 1]);  % Var hyperparam drawn from MVGamma(shape=5/2, scale=5/2)
Sig_0 = diag(Sig_0);                        % Assume indep components. (pcov)x(pcov) matrix of variances. 

% Prior for model params, drawn from MVN(mu0, Sig0)
probit_params.xi = mvnrnd(mu_0, Sig_0);  

% Probit model latent variable z_i 
probit_params.z_i = zeros(data_vars.n, 1);   % Initialize latent probit variable, z_i
% Factor variable coding design matrix with demographic covariates and initial latent classes
Q = zeros(data_vars.n, S * k_fixed);
for s = 1:S
    for k = 1:k_fixed
        Q(:, ((s-1) * k_fixed) + k) = q_dem(:, s) .* OFMM_params.x_ci(:, k);
    end
end
lin_pred = Q * transpose(probit_params.xi);  % Linear predictor, Q*xi. Mean of truncated normal dist
% For cases, z_i ~ TruncNormal(mean=Q*xi, var=1, low=0, high=Inf)
probit_params.z_i(data_vars.y == 1) = truncnormrnd(1, lin_pred(data_vars.y == 1), 1, 0, Inf); 
% For controls, z_i ~ TruncNormal(mean=Q*xi, var=1, low=-Inf, high=0)
probit_params.z_i(data_vars.y == 0) = truncnormrnd(1, lin_pred(data_vars.y == 0), 1, -Inf, 0); 
if sum(probit_params.z_i == Inf) > 0                  % Control extremes
    probit_params.z_i(probit_params.z_i == Inf) = norminv(1 - 1e-10);
end 
if sum(probit_params.z_i == -Inf) > 0
    probit_params.z_i(probit_params.z_i == -Inf) = norminv(1e-10);
end 

% Individual probit log-likelihood with class assignments and latent variable formulation    
probit_params.probit_lik = normpdf(probit_params.z_i, lin_pred, 1) .* ((data_vars.y == 1).* (probit_params.z_i > 0) + (data_vars.y == 0) .* (probit_params.z_i <= 0));

% Probit likelihood component of posterior class membership P(c_i|-) with assumed classes
probit_params.indiv_lik_probit_class = zeros(data_vars.n, k_max); % Initialize matrix of indiv probit likelihood assuming each class       
for k = 1:k_max
    % Factor variable coding temporarily assign all indivs to class k
    Q_temp = zeros(data_vars.n, S * k_max);  % Initialize temp design matrix for class membership
    for s = 1:S
        Q_temp(:, ((s-1) * k_max) + k) = q_dem(:, s);         % Temporarily assign all indivs to class k
    end
    lin_pred_temp = Q_temp * transpose(probit_params.xi); % Linear predictor, Q*xi
    % Update indiv probit likelihood assuming class k
    probit_params.indiv_lik_probit_class(:, k) = normpdf(probit_params.z_i,lin_pred_temp,1).*((data_vars.y==1).*(probit_params.z_i>0)+(data_vars.y==0).*(probit_params.z_i<=0));
end  


% Run MCMC algorithm using fixed number of classes
n_runs = 25000;  % Number of MCMC iterations
burn = 15000;    % Burn-in period
thin = 5;        % Thinning factor
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
    % Update posterior class membership probs, pi
    for k = 1:k_max  % For each latent class
        % Get updated num indivs assigned to class k
        OFMM_params.n_ci(k) = sum(OFMM_params.c_i == k);  
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
        prob_ci_numer(:, k) = OFMM_params.pi(k) * prod(categ_lik, 2) .* probit_params.indiv_lik_probit_class(:, k);  
    end
    % Matrix of updated posterior P(c_i=k|-) for each indiv and class
    prob_ci = bsxfun(@times, prob_ci_numer, 1 ./ (sum(prob_ci_numer, 2)));
    % Adjust precision for mnrnd function
    for i = 1:data_vars.n
       prob_i = prob_ci(i, :);
       if any(prob_i < 1e-15)
           prob_i = prob_i + 1e-15;
           prob_ci(i, :) = prob_i ./ sum(prob_i);
       end    
    end       
    OFMM_params.x_ci = mnrnd(1, prob_ci);                   % Matrix of updated c_i drawn from Mult(1, prob_ci). Each row is draw for an indiv
    [row, col] = find(OFMM_params.x_ci == 1);               % Row and col indices of each nonzero element of x_ci; col is class assignment for each indiv
    sorted = sortrows([row col], 1);            % Sort indices in ascending row index order, to match subject index
    OFMM_params.c_i = sorted(:, 2);             % Vector of updated class assignment, c_i, for each individual
      
    
    % Update posterior item response probabilities, theta   
    n_theta = zeros(data_vars.p, data_vars.d_max); % Initialize matrix of num indivs with each item response level
    for k = 1:k_max  
        c_i_mat = repmat(OFMM_params.c_i, [1, data_vars.p]);  % Rep c_i's p times to form matrix of classes for each indiv and item
        class_wt = (c_i_mat == k);                            % Matrix indicating indivs assigned to class k. All other rows set to 0
        for r = 1:data_vars.d_max                             % For each consumption level r
            
            % DEBUGGING ERROR
            if any(size((data_vars.food == r)) ~= size(class_wt))
                disp("FOOD");
                disp(size((data_vars.food == r)));
                disp("CLASS_WT");
                disp(size(class_wt));
                disp("SORTED");
                disp(size(sorted));
                disp(sorted);
            end
            
            % Vector of num indivs with level r and class k, for each item (i.e., sum I(x_ij=r, c_i=k) over indivs)
            n_theta(:, r) = sum((data_vars.food == r) .* class_wt)'; 
        end
        for j = 1:data_vars.p                                     % For each food item
            dj = data_vars.d(j);                                  % Num response levels for item j
            eta_post = eta(1:dj) + n_theta(j, 1:dj);              % Posterior hyperparam is eta[r] + n_theta(j,r)
            OFMM_params.theta(j, k, 1:dj) = drchrnd(eta_post, 1); % Draw from posterior dist: Dir(eta), for every item and class
        end
    end        
    
    
    % Update latent probit variable, z_i
    % Factor variable coding design matrix with demographic covariates and initial latent classes
    Q = zeros(data_vars.n, S * k_max);
    for s = 1:S
        for k = 1:k_max
            Q(:, ((s-1) * k_max) + k) = q_dem(:, s) .* OFMM_params.x_ci(:, k);
        end
    end
    lin_pred = Q * transpose(probit_params.xi);  % Linear predictor, Q*xi. Mean of truncated normal dist
    % For cases, z_i ~ TruncNormal(mean=Q*xi, var=1, low=0, high=Inf)
    probit_params.z_i(data_vars.y == 1) = truncnormrnd(1, lin_pred(data_vars.y == 1), 1, 0, Inf); 
    % For controls, z_i ~ TruncNormal(mean=Q*xi, var=1, low=-Inf, high=0)
    probit_params.z_i(data_vars.y == 0) = truncnormrnd(1, lin_pred(data_vars.y == 0), 1, -Inf, 0); 
    if sum(probit_params.z_i == Inf) > 0                       % Control extremes
        probit_params.z_i(probit_params.z_i == Inf) = norminv(1 - 1e-10);
    end 
    if sum(probit_params.z_i == -Inf) > 0
        probit_params.z_i(probit_params.z_i == -Inf) = norminv(1e-10);
    end    
    probit_params.probit_lik = normpdf(probit_params.z_i, lin_pred, 1) .* ((data_vars.y == 1).* (probit_params.z_i > 0) + (data_vars.y == 0) .* (probit_params.z_i <= 0));

    
    % Update individual log-likelihood with class assignments and latent variable formulation       
    item_idx = repmat(1:data_vars.p, data_vars.n, 1);     % nxp matrix of item id's. Replicates 1:p row n times
    item_idx = item_idx(:);                               % Concat by col into (np)x1 vector. 1(x n),2(x n),...,p(x n)
    class_idx = repmat(OFMM_params.c_i, 1, data_vars.p);  % Replicate col vector c_i p times to form nxp matrix of class assigns
    class_idx = class_idx(:);                             % Concat by col into (np)x1 vector 
    x_idx = data_vars.food(:);                            % Concat responses by col in (np)x1 vector 
    % lin_idx: (npx1) vector of linear indices indicating value of pxKxd theta matrix, for each item and indiv
    lin_idx = sub2ind([data_vars.p, k_max, data_vars.d_max], item_idx, class_idx, x_idx);
    % nxp matrix of theta values for each indiv and item \prod_{r=1}^d \theta_{jr|k}^{I(x_ij=r, c_i=k)}
    theta_indiv = reshape(OFMM_params.theta(lin_idx), [data_vars.n, data_vars.p]);
    % Indiv complete log-likelihood
    indiv_loglik = log(OFMM_params.pi(OFMM_params.c_i) * prod(theta_indiv, 2) * probit_params.probit_lik);
    
    
    % Update posterior probit model coefficients, xi   
    Sig_post = inv(Sig_0) + (transpose(Q) * Q);           % Precision of posterior Normal dist for xi
    mu_right = (Sig_0 \ mu_0) + (transpose(Q) * probit_params.z_i);       % Right part of mean of posterior Normal dist for xi
    % mu_right = (Sig_0\probit_params.xi0') + (transpose(Q)*z_i);  CHANGE THIS
    mu_post = Sig_post \ mu_right;                      % Mean of posterior Normal dist for xi
    probit_params.xi = mvnrnd(mu_post, inv(Sig_post));  % Draw from posterior dist: N(mu_post, Sig_post^{-1})    
    
    
    % Update probit likelihood component of posterior class membership P(c_i|-)
    for k = 1:k_max
        % Factor variable coding temporarily assign all indivs to class k
        Q_temp = zeros(data_vars.n, S * k_max);  % Initialize temp design matrix for class membership
        for s = 1:S  % (s-1) * k_max is the shift amount
            Q_temp(:, ((s-1) * k_max) + k) = q_dem(:, s);         % Temporarily assign all indivs to class k
        end
        lin_pred_temp = Q_temp * transpose(probit_params.xi); % Linear predictor, Q*xi       
        % Update indiv probit likelihood assuming class k
        probit_params.indiv_lik_probit_class(:, k) = normpdf(probit_params.z_i,lin_pred_temp,1).*((data_vars.y==1).*(probit_params.z_i>0)+(data_vars.y==0).*(probit_params.z_i<=0));
    end    
    
    
    % Store posterior values if the iteration is selected based on thinning  
    if mod(iter, thin) == 0
        MCMC_out.pi(iter / thin, :) = OFMM_params.pi;
        MCMC_out.c_i(iter / thin, :) = OFMM_params.c_i;
        MCMC_out.theta(iter / thin, :, :, :) = OFMM_params.theta; 
        MCMC_out.loglik(iter / thin) = sum(indiv_loglik); 
        MCMC_out.z_i(iter / thin, :) = probit_params.z_i;
        MCMC_out.xi(iter / thin, :) = probit_params.xi; 
    end
    
    % Relabel classes every 10 iterations to encourage mixing
    if mod(iter, 10) == 0
        new_order = randperm(k_max);                               % New labels for latent classes
        new_ci = OFMM_params.c_i;
        for k = 1:k_max
            new_ci(OFMM_params.c_i == k) = new_order(k);           % Class assignments with new labels
        end
        OFMM_params.c_i = new_ci;                               % Relabel class assignments
        OFMM_params.theta = OFMM_params.theta(:, new_order, :); % Relabel item-response probabilities
        % Relabel indiv probit likelihood
        probit_params.indiv_lik_probit_class = probit_params.indiv_lik_probit_class(:, new_order);  
        OFMM_params.pi = OFMM_params.pi(new_order);                % Relabel class probabilities
    end    
end
MCMC_out.runtime = toc;  % Stop timer

% Discard burn-in
MCMC_out.pi = MCMC_out.pi((burn / thin) + 1:end, :);
MCMC_out.theta = MCMC_out.theta((burn / thin) + 1:end, :, :, :);
MCMC_out.c_i = MCMC_out.c_i((burn / thin) + 1:end, :);
MCMC_out.xi = MCMC_out.xi((burn / thin) + 1:end, :);  
MCMC_out.z_i = MCMC_out.z_i((burn / thin) + 1:end, :);
MCMC_out.loglik = MCMC_out.loglik((burn / thin) + 1:end); 


% Post-processing to recalibrate labels and remove extraneous empty classes
% Remove extraneous empty classes and relabel class assignments  
m = size(MCMC_out.pi, 1);                                 % Num stored output iterations
post_MCMC_out.k_med = median(sum(MCMC_out.pi > 0.05, 2)); % Median num classes w/ more than >5% indivs, over all iterations
% Hierarchical clustering dendrogram tree using pairwise distance of num times two indivs were in same class together
post_MCMC_out.tree = linkage(single(MCMC_out.c_i'), 'complete', 'hamming');
% Vector of new class assignments for each indiv
red_ci = cluster(post_MCMC_out.tree, 'maxclust', post_MCMC_out.k_med); 
relabel_ci = zeros(m, post_MCMC_out.k_med);               % Initialize matrix of relabeled class assigns after removing empty classes
for k = 1:post_MCMC_out.k_med                             % For each non-empty class
    % For each MC run, among indivs assigned to new class k, which original class was most common 
    relabel_ci(:, k) = mode(MCMC_out.c_i(:, red_ci == k), 2); 
end

% Reorder parameter estimates to match the correct new class assignment
post_MCMC_out.pi = zeros(m, post_MCMC_out.k_med); % Initialize parameter estimates corresponding to the non-empty classes
post_MCMC_out.theta = zeros(m, data_vars.p, post_MCMC_out.k_med, data_vars.d_max);
% Factor variable initialize xi
post_MCMC_out.xi = zeros(m, S * post_MCMC_out.k_med);
k_old = size(MCMC_out.pi, 2);
for iter = 1:m                                    % For each stored MC run
    iter_order = relabel_ci(iter, :);             % Vector of new class assignments (mode orig classes after removing empty classes)
    pi_temp = MCMC_out.pi(iter, iter_order);      % Update class probs corresp to new class assignments
    post_MCMC_out.pi(iter, :) = pi_temp / sum(pi_temp);  % Normalize to sum to 1
    % Update item-response probs corresp to new class assignments
    post_MCMC_out.theta(iter, :, :, :) = MCMC_out.theta(iter, :, iter_order, :);   
    % Factor variable coding update probit model coefs corresp to new class assignments
    for s = 1:S
        post_MCMC_out.xi(iter, (s-1)*post_MCMC_out.k_med + (1:post_MCMC_out.k_med)) = MCMC_out.xi(iter, (s-1)*k_old + iter_order);
    end
%             % Update probit model coefs corresp to new class assignments. First few coefs corresp to dem vars and are not changed
%             post_MCMC_out.xi(iter, :) = [MCMC_out.xi(iter, 1:S) MCMC_out.xi(iter, S + iter_order)]; 
end
post_MCMC_out.loglik = MCMC_out.loglik;  % Obtain log-lik for each stored MC run


%Obtain posterior estimates, reduce number of classes, analyze results, and save output
% Identify unique classes   
% Array of posterior median item-response probs
theta_med_temp = reshape(median(post_MCMC_out.theta), [data_vars.p, post_MCMC_out.k_med, data_vars.d_max]);  
% Matrix of most likely consump level based on item-response probs, for each item and class across MCMC runs
[~, ind0] = max(theta_med_temp, [], 3); 
t_ind0 = transpose(ind0);
[u_mat, ia, ~] = unique(t_ind0, 'rows');  % Vector of class indices corresponding to unique classes
analysis.k_red = length(ia);          % Number of unique classes

% Posterior MCMC chains for pi after combining duplicated classes, re-normalized
analysis.pi_red = post_MCMC_out.pi(:, ia);    % Restrict to unique classes
has_dupes = size(u_mat, 1) < size(t_ind0, 1); % Check for duplicates
if has_dupes
    ind_dup = setdiff(1:size(t_ind0,1), ia);  % Find indices of duplicated classes
    for i = 1:length(ind_dup)                    % Iterate over dupe rows
        for j = 1:length(ia)                     % Iterate over non-dupe rows
            i_dup = ind_dup(i);               
            i_orig = ia(j);
            if t_ind0(i_dup, :) == t_ind0(i_orig, :)
                % Combine pi for duplicated patterns
                analysis.pi_red(:, i_orig) = analysis.pi_red(:, i_orig) + post_MCMC_out.pi(:, i_dup);
                %disp(median(analysis.pi_red));
            end
        end
    end
end  
analysis.pi_red = analysis.pi_red ./ sum(analysis.pi_red, 2);

% Posterior MCMC parameter chains for the unique classes, re-normalized
analysis.theta_red = post_MCMC_out.theta(:, :, ia, :);
analysis.theta_red = analysis.theta_red ./ sum(analysis.theta_red, 4);
% Factor variable version
analysis.xi_red = zeros(size(analysis.pi_red, 1), S * analysis.k_red);
for s = 1:S
    analysis.xi_red(:, ((s-1)*analysis.k_red) + (1:analysis.k_red)) = post_MCMC_out.xi(:, (s-1)*post_MCMC_out.k_med + ia);
end
%         analysis.xi_red = post_MCMC_out.xi(:, [1:S (S+ia')]);    

% Posterior median estimates for the unique classes
analysis.pi_med = median(analysis.pi_red);  % Class probs 
analysis.pi_med = analysis.pi_med / sum(analysis.pi_med);
% Item-response probs 
analysis.theta_med = reshape(median(analysis.theta_red), [data_vars.p, analysis.k_red, data_vars.d_max]);  
analysis.theta_med = analysis.theta_med ./ sum(analysis.theta_med, 3);
analysis.xi_med = median(analysis.xi_red);  % Probit coefs      

% Update class membership probit component using unique classes and updated posterior median estimates for xi
indiv_lik_probit_class = zeros(data_vars.n, analysis.k_red); % Initialize indiv probit lik matrix for each indiv and class
analysis.z_i = MCMC_out.z_i(end, :)';                        % Obtain last values of probit latent variable
for k = 1:analysis.k_red
    % Factor variable coding temporarily assign all indivs to class k
    Q_temp = zeros(data_vars.n, S * analysis.k_red);  % Initialize temp design matrix for class membership
    for s = 1:S
        Q_temp(:, ((s-1) * analysis.k_red) + k) = q_dem(:, s);         % Temporarily assign all indivs to class k
    end
    lin_pred_temp = Q_temp * transpose(analysis.xi_med); % Linear predictor, Q*xi
    % Update indiv probit likelihood assuming class k
    indiv_lik_probit_class(:, k) = normpdf(analysis.z_i,lin_pred_temp,1).*((data_vars.y==1).*(analysis.z_i>0)+(data_vars.y==0).*(analysis.z_i<=0));
end  

% Update class membership, c_i, using unique classes and updated posterior median estimates for pi and theta
prob_ci_numer = zeros(data_vars.n, analysis.k_red);  % Initialize indiv OFMM lik matrix for each indiv and class
for k = 1:analysis.k_red
    % Matrix of unique item-level probs for class k
    theta_k = reshape(analysis.theta_med(:, k, :), data_vars.p, data_vars.d_max); 
    % Matrix of categ likelihood \prod_{r=1}^d \theta_{jr|k}^{I(x_ij=r)} for each indiv and item 
    categ_lik = reshape(theta_k(data_vars.lin_idx), [data_vars.n, data_vars.p]); 
    % Update numerator for posterior P(c_i=k|-) (indiv OFMM likelihood)
    prob_ci_numer(:, k) = analysis.pi_med(k) * prod(categ_lik, 2) .* indiv_lik_probit_class(:, k);   
end
% Matrix of updated posterior P(c_i=k|-) for each indiv and class
prob_ci = bsxfun(@times, prob_ci_numer, 1 ./ (sum(prob_ci_numer, 2)));
% Max posterior class assign probs and class assigns for each indiv using posterior median estimates
[analysis.max_prob_ci, analysis.c_i] = max(prob_ci, [], 2);

% Obtain probit model mean using median parameter estimates
% Factor variable coding design matrix with updated class memberships using median estimates 
Q_med = zeros(data_vars.n, S * analysis.k_red);
x_ci = dummyvar(analysis.c_i);
for s = 1:S
    for k = 1:analysis.k_red
        Q_med(:, ((s-1) * analysis.k_red) + k) = q_dem(:, s) .* x_ci(:, k);
    end
end
analysis.Phi_med = normcdf(Q_med * transpose(analysis.xi_med)); % Linear predictor, Q*xi. Probit model mean using median estimates
analysis.Phi_med(analysis.Phi_med == 0) = 1e-10;  % Adjust extremes
analysis.Phi_med(analysis.Phi_med == 1) = 1 - 1e-10;    

% Update individual log-likelihood using unique classes and updated posterior median estimates 
% Categorical component
item_idx = repmat(1:data_vars.p, data_vars.n, 1);     % nxp matrix of item id's. Replicates 1:p row n times
item_idx = item_idx(:);                               % Concat by col into (np)x1 vector. 1(x n),2(x n),...,p(x n)
class_idx = repmat(analysis.c_i, 1, data_vars.p);     % Replicate col vector c_i p times to form nxp matrix of class assigns
class_idx = class_idx(:);                             % Concat by col into (np)x1 vector 
x_idx = data_vars.food(:);                            % Concat responses by col in (np)x1 vector 
% lin_idx: (npx1) vector of linear indices indicating value of pxKxd theta matrix, for each item and indiv
lin_idx = sub2ind([data_vars.p, analysis.k_red, data_vars.d_max], item_idx, class_idx, x_idx);
% nxp matrix of theta values for each indiv and item \prod_{r=1}^d \theta_{jr|k}^{I(x_ij=r, c_i=k)}
theta_indiv = reshape(analysis.theta_med(lin_idx), [data_vars.n, data_vars.p]);
% Probit component 
probit_lik = normpdf(analysis.z_i, analysis.Phi_med, 1) .* ((data_vars.y == 1).* (analysis.z_i > 0) + (data_vars.y == 0) .* (analysis.z_i <= 0));
% Indiv complete log-likelihood
analysis.loglik_med = log(analysis.pi_med(analysis.c_i) * prod(theta_indiv, 2) * probit_lik);

% DIC is a metric for model goodness of fit. It consists of two terms:
    % 1) median(loglik): calculate log-lik at each MCMC iter then get the median
    % 2) loglik_med: get posterior estimates of each param, then use these to calculate the log-lik
% DIC-6 is an adaptation that penalizes overfitting
analysis.dic6 = -6 * median(post_MCMC_out.loglik) + 4 * sum(analysis.loglik_med); 

% AEBIC is a metric for model goodness of fit that is asymptotically
% consistent and invariant to reparameterization
num_params = numel(analysis.pi_med) + numel(analysis.theta_med) + numel(analysis.xi_med);
gamma = 1;
t1 = -2 * data_vars.n * mean(post_MCMC_out.loglik ./ data_vars.n); % -2*n*E[1/n*sum(indiv_log_lik)]
t2 = num_params * log(data_vars.n);                                % |S|*log(n)
t3 = 2 * gamma * num_params * log(p_cov);                          % 2*gamma*|S|*log(p)
analysis.aebic = t1 + t2 + t3;    
    

%% CHECK CONVERGENCE
figure; %check mixing of pi parameter
plot(MCMC_out.pi)

figure; %check ordered pi
plot(post_MCMC_out.pi) 

figure; % check ordered theta
plot(post_MCMC_out.theta(:,1,1,1))
hold on
plot(post_MCMC_out.theta(:,1,1,2))
hold off

figure; % check ordered xi
plot(post_MCMC_out.xi)

figure; % check ordered theta
plot(post_MCMC_out.theta(:,1,2,1))
hold on
plot(post_MCMC_out.theta(:,1,2,2))
hold off

immse(samp_data.true_Phi, analysis.Phi_med)
sum(abs(sort(samp_data.true_Phi) - sort(analysis.Phi_med)))

%% CONVERTING FACTOR VARIABLE TO REFERENCE CELL
if samp_data.true_Ci(1) ~= analysis.c_i(1)
    beta_1 = analysis.xi_med(3);
    beta_2 = analysis.xi_med(4) - analysis.xi_med(3);
    beta_3 = analysis.xi_med(1) - analysis.xi_med(3);
    beta_4 = analysis.xi_med(2) - beta_1 - beta_2 - beta_3;  
else
    beta_1 = analysis.xi_med(4);
    beta_2 = analysis.xi_med(3) - analysis.xi_med(4);
    beta_3 = analysis.xi_med(2) - analysis.xi_med(4);
    beta_4 = analysis.xi_med(1) - beta_1 - beta_2 - beta_3;
end
beta = [beta_1 beta_2 beta_3 beta_4]
% Find individuals with the four patterns
Q1 = find((samp_data.true_Si == 2) & (samp_data.true_Ci == 2));
Q2 = find((samp_data.true_Si == 1) & (samp_data.true_Ci == 2));
Q3 = find((samp_data.true_Si == 2) & (samp_data.true_Ci == 1));
Q4 = find((samp_data.true_Si == 1) & (samp_data.true_Ci == 1));
Q_indivs = [Q1(1); Q2(1); Q3(1); Q4(1)];
pred_Phi = analysis.Phi_med(Q_indivs)
%pred_Phi = normcdf(Q_test_ref * transpose(beta));
% Factor variable
Q_test_ref = [1 0 0 0;   % S=1,C=1
              0 1 0 0;   % S=1,C=2
              0 0 1 0;   % S=2,C=1
              0 0 0 1];  % S=2,C=2
%     Q_test_ref = [1 0 0 0;   % S=2,C=2
%                   1 1 0 0;   % S=1,C=2
%                   1 0 1 0;   % S=2,C=1
%                   1 1 1 1];  % S=1,C=1
actual_Phi = normcdf(Q_test_ref * transpose(samp_data.true_xi));
sum(abs(sort(actual_Phi) - sort(pred_Phi)))
normcdf(Q_test_ref * transpose(beta))


% Find individuals with incorrect class assignment
sim_data = samp_data;
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
sens = mean(sens_class) % Mean sensitivity over all true classes
spec = mean(spec_class) % Mean specificity over all true classes

% Find incorrectly classified individuals
if samp_data.true_Ci(1) ~= analysis.c_i(1)
    temp_ci = ones(size(analysis.c_i));
    temp_ci(analysis.c_i == 1) = 2;
    mismatch = find(temp_ci ~= samp_data.true_Ci);
else
    mismatch = find(analysis.c_i ~= samp_data.true_Ci);
end    
length(mismatch)



% Find differences using observed Phi values
true_Phi_values = zeros(length(samp_data.true_xi), 1);
% P(Y=1|S=1,C=1)
subset = samp_data.Y_data(samp_data.true_Si == 1 & samp_data.true_Ci == 1);
true_Phi_values(1) = sum(subset) / length(subset);
% P(Y=1|S=1,C=2)
subset = samp_data.Y_data(samp_data.true_Si == 1 & samp_data.true_Ci == 2);
true_Phi_values(2) = sum(subset) / length(subset);
% P(Y=1|S=2,C=1)
subset = samp_data.Y_data(samp_data.true_Si == 2 & samp_data.true_Ci == 1);
true_Phi_values(3) = sum(subset) / length(subset);
% P(Y=1|S=2,C=2)
subset = samp_data.Y_data(samp_data.true_Si == 2 & samp_data.true_Ci == 2);
true_Phi_values(4) = sum(subset) / length(subset);
true_Phi_values
sum(abs(sort(true_Phi_values) - sort(pred_Phi)))

% sum abs dev and sum squared distance for xi's
sum(abs(sort(samp_data.true_xi) - sort(analysis.xi_med)))
sum((sort(samp_data.true_xi) - sort(analysis.xi_med)).^2) % take this over iterations to get MSE

% si = dummyvar([1 1 2 2 2]);
% ci = dummyvar([1 2 3 3 2]);
% S = size(si, 2);
% K = size(ci, 2);
% Q = zeros(size(si, 1), S*K);
% shift = 0;
% for s = 1:S
%     for k = 1:K
%         Q(:, shift + k) = si(:, s) .* ci(:, k);
%     end
%     shift = shift + K;
% end
% Q
% Q = zeros(size(si, 1), S*K);
% shift = 0;
% k = 1;
% for s = 1:S
%     Q(:, shift + k) = si(:, s);
%     shift = shift + K;
% end
% Q
% xi_orig = [1 2 3 4 5 6 7 8];  % orig 2 subpop, 4 classes
% order = [4 2];
% k_old = 4;
% k_new = 2;
% S = 2;
% xi_new = zeros(1, 2 * k_new);
% shift = 0;
% for s = 1:S
%     xi_new(((s-1)*k_new) + (1:k_new)) = xi_orig(((s-1)*k_old) + order);
% end
% xi_new