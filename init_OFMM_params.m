% init_OFMM_params initializes priors and variables for the OFMM model given hyperparameters
% Inputs: 
%   data_vars: output from wtd_get_data_vars or get_data_vars function
%   k_max: max number of latent classes
%   alpha: vector of hyperparams for class membership probs, pi; 1x(k_max). Default is alpha=1/50
%   eta: vector of hyperparams for item response probs, theta; 1x(d_max). Default is eta=1
% Outputs: OFMM_params structural array with the following fields:
%   pi: vector prior for class membership probs; 1xp. Dist: Dir(alpha)
%   theta: vector prior for item response probabilities; px(k_max)x(d_max). Dist: Dir(eta) per class and item
%   c_i: vector of initial class assignments (randomly generated); nx1
%   n_ci: vector of num indivs initially assigned to each cluster; 1x(k_max)
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
    OFMM_params.n_ci=sum(x_ci);               % Vector of num indivs assigned to each cluster
end