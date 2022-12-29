% post_process_unsup_latent conducts post-processing to remove extraneous 
% empty classes and relabel class assignments for the wOFMM model
% Inputs: 
%   MCMC_out: output initialized in wtd_run_MCMC_latent or run_MCMC_latent function
%   data_vars: output from wtd_get_data_vars_latent or get_data_vars_latent function
% Outputs: returns post_MCMC_out structural array with the following fields:
%   k_med: Median num classes w/ more than >5% indivs, calculated over MCMC iterations
%   tree: Hierarchical clustering dendrogram tree using distance criterion to find num classes
%   pi: Matrix of reordered class membership probs; (n_runs/thin)x(k_med)
%   theta: Matrix reordered item-response probs; (n_runs/thin)xpx(k_med)x(d_max)
function post_MCMC_out = post_process_unsup_latent(MCMC_out, data_vars) 
    % Remove extraneous empty classes and relabel class assignments  
    m = size(MCMC_out.pi, 1);                                 % Num stored output iterations
    post_MCMC_out.k_med = round(median(sum(MCMC_out.pi > 0.05, 2))); % Median num classes w/ more than >5% indivs, over all iterations
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
    for iter = 1:m                                    % For each stored MC run
        iter_order = relabel_ci(iter, :);             % Vector of new class assignments (mode orig classes after removing empty classes)
        pi_temp = MCMC_out.pi(iter, iter_order);      % Update class probs corresp to new class assignments
        post_MCMC_out.pi(iter, :) = pi_temp / sum(pi_temp);  % Normalize to sum to 1
        % Update item-response probs corresp to new class assignments
        post_MCMC_out.theta(iter, :, :, :) = MCMC_out.theta(iter, :, iter_order, :);
    end
end