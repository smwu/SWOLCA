% init_OFMM_params_latent initializes priors and variables for the OFMM model given hyperparameters
% Inputs: 
%   data_vars: output from wtd_get_data_vars
%   k_max: max number of latent classes
%   alpha: vector of hyperparams for class membership probs, pi; 1x(k_max). Default is alpha=1/50
%   eta: vector of hyperparams for item response probs, theta; 1x(d_max). Default is eta=1
% Outputs: OFMM_params structural array with the following fields:
%   pi: vector prior for class membership probs; 1xp. Dist: Dir(alpha)
%   theta: vector prior for item response probabilities; px(k_max)x(d_max). Dist: Dir(eta) per class and item
%   c_i: vector of initial class assignments (randomly generated); nx1
%   n_ci: vector of num indivs initially assigned to each cluster; 1x(k_max)
%   x_ci: matrix of class assignments for all indivs; nx(k_max)
function OFMM_params = wtd_init_OFMM_params_latent(data_vars, k_max, alpha, eta)
    % Prior for class membership probs
    OFMM_params.pi = drchrnd(alpha, 1); 
    
    % Random initialization of class assignments
    OFMM_params.x_ci = mnrnd(1, OFMM_params.pi, data_vars.n); % Matrix of c_i drawn from Mult(1, pi). Each row is draw for an indiv
    [row, col] = find(OFMM_params.x_ci);                      % Row and col indices of nonzero elements of x_ci; col is class assign for each indiv
    sorted = sortrows([row col], 1);              % Sort indices in ascending row index order, to match subject index
    OFMM_params.c_i = sorted(:, 2);               % Vector of class assignment, c_i, for each individual
    OFMM_params.n_ci = zeros(1, k_max);           % Initialize vector of weighted num indivs assigned to each class
    for k = 1:k_max                               % For each latent class
        % Get weighted num indivs assigned to class k
        OFMM_params.n_ci(k) = sum(data_vars.wt_kappa(OFMM_params.c_i == k));  
    end
    
    % Prior for item response probabilities
    OFMM_params.theta = zeros(data_vars.p, k_max, data_vars.d_max);  % Initialize array of item response probs
    for k = 1:k_max                % For each class
        for j = 1:data_vars.p      % For each food item
            d_j = data_vars.d(j);  % Max number of levels for food item j
            OFMM_params.theta(j, k, 1:d_j) = drchrnd(eta(1:d_j), 1); % For each class and item, draw from Dir(eta)
        end
    end     
end