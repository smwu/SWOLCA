

% Each MCMC iteration

%% Get data variables
% Outputs: data_vars structural array with the following fields:
%   food: matrix of food intake data
%   n: number of individuals
%   p: number of food items
%   d_max: max number of consumption levels over all items
%   d: vector of max number of levels for each food item; px1
%   y: vector of outcomes; nx1
%   wt: vector of survey weights; nx1
%   wt_kappa: vector of normalized weights; nx1
%   wt_kappa_mat: matrix of normalized weights, replicated across items; nxp
%   lin_idx: vector of linear indices for unique item-response combos; (n*p)x1
data_vars = wtd_get_data_vars_latent(samp_data);

% % Normalized standardized weights that sum to 1
% norm_wt_kappa = data_vars.wt_kappa ./ sum(data_vars.wt_kappa);
% combined_data = [samp_data.X_data, samp_data.Y_data, samp_data.sample_wt];

test.X_data = zeros(5, 1);
test.Y_data = [1; 2; 3; 4; 5];
weights = [10; 20; 30; 40; 50];
N = sum(weights);
n = length(test.Y_data);

tic;
M = 10000;
y_bar = zeros(M, 1);
for m = 1:M
    synth_pop = wt_polya_post(test, weights, n, N, m);
    y_bar(m) = sum(synth_pop.Y_data) / N;
end
disp(mean(y_bar));
toc
% Truth: 3.6666
% M=100: 3.6195 (0.24 secs)
% M=1000: 3.6444 (0.29 secs)
% M=10000: 3.6455 (1.37 secs)

