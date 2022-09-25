% wtd_get_data_vars takes in sample data and outputs relevant variables for
% the wsOFMM model
% Inputs: samp_data structural array with at least the following columns:
%   X_data: food intake data as a matrix
%   Y_data: outcome data as a vector
%   sample_wt: survey weights as a vector
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
function data_vars = wtd_get_data_vars(samp_data)
    data_vars.food = samp_data.X_data;
    [data_vars.n, data_vars.p] = size(data_vars.food);
    data_vars.d_max = max(data_vars.food(:));  % Max number of levels over all items
    data_vars.d = max(data_vars.food);         % Max number of levels for each food item. 
    data_vars.y = samp_data.Y_data;            % Outcome
    data_vars.wt = samp_data.sample_wt;        % Survey weights    
    
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
end