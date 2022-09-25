% wt_polya_post uses a weighted Polya posterior to generate a synthetic
% population, given the sample data, weights, sample and population sizes,
% and a random seed.
% Inputs:
%   X_data: food intake data as a matrix; nxp
%   weights: Vector of weights; nx1
%   n: Sample size
%   N: Population size
%   seed: random number generator seed
%   varargin: additional variable for Y data vector; nx1
% Outputs:
%   synth_pop: Structural array with the following elements:
%       X_data: synthetic population data imputed for food intake data
%       Y_data: synthetic population data imputed for outcome data
function synth_pop = wt_polya_post(X_data, weights, n, N, seed, varargin)

    rng(seed, 'twister');
    w = cumsum(weights);
    
    imp_index = zeros(N - n, 1); % Initialize vector of imputed indices
    for i = 1:(N-n)              % For each non-sampled unit
        % Generate a weighted Polya draw
        a = w(n) * rand;         % Random number from 0 to sum(weights)
        j = 1;                   % j is first index with a <= cumsum
        while a > w(j)
            j = j + 1;
        end
        % To represent adding another ball of the same color, increment all
        % cumsum values starting from the j-th value
        for k = j:n                
            w(k) = w(k) + 1;
        end
        % Index indicates i-th imputed unit has same value as j-th sample unit
        imp_index(i) = j; 
    end
    
    % Synthetic population combines the sample data with imputed data
    synth_pop.X_data = [X_data; X_data(imp_index, :)];
    % Combine Y data if provided as an input
    if (nargin > 5)
        Y_data = varargin{1};
        synth_pop.Y_data = [Y_data; Y_data(imp_index, :)];
    end

end

%% Old Code
%     l_boot = zeros(n, 1); % Num selections for each obs among previous draws
%     d = 0;                % Index of draw
%     Nnn = (N - n) / n;
% 
%     % Weighted Polya sampling is equivalent to selecting an observation with 
%     % probability weight (w_i^* + l_{i,d-1}) / (n+d-1) at each draw d, where  
%     %   w_i^* = n/(N-n)*(w_i-1) is the normalized indiv wt, normalized to 
%     % sum to n and with certainty units excluded. 
%     %   w_i = original wt (replicate weight if bootstrap sample)
%     %   l_{i,d-1} = number of selections of the obs among previous d-1 draws
%     %   d = index of draw
%     % Equiv to Cohen (1997) formulation
%     % (w_i - 1 + l_{i,d-1}(N-n)/n) / (N - n + (d-1)(N-n)/n)) 
%     wt_cohen = (weights - 1 + l_boot * Nnn) ./ (N - data_vars.n + (d - 1) * Nnn); 
% 
%     % Initialize matrix of consumption exposure data for non-sampled indivs
%     impute.X_data = zeros(N - n, data_vars.p); 
%     % Initialize vector of CVD response data for non-sampled indivs
%     impute.Y_data = zeros(N - data_vars.n, 1);
%     for m_ind = 1:m  % For each non-sampled individual
%         % Bootstrap sample 1 item with replacement from sample data, according to Polya weights
%         % y_select = sampled values. Sampled value is assigned to kth value of y_select
%         % idx = sampled index
%         imp_index = datasample(1:data_vars.n, 1, 'Replace', true, 'Weights', wt_cohen);
%         impute.X_data(m_ind, :) = samp_data.X_data(imp_index, :);
%         impute.Y_data(m_ind) = samp_data.Y_data(imp_index);
% 
%         % Increment variables
%         l_boot(m_ind) = l_boot(m_ind) + 1;
%         d = d + 1;
%     end
%     
%     % Synthetic population combines the sample data with the bootstrapped non-sampled data
%     synth_pop.X_data = [samp_data.X_data; impute.X_data];
%     synth_pop.Y_data = [samp_data.Y_data; impute.Y_data];
%     % Form PRS by taking a SRS from the population
%     PRS_index = datasample(1:N, data_vars.n, 'Replace', false);