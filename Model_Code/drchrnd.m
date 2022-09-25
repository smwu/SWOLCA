% drchrnd generates a random sample from a Dirichlet distribution
% Inputs:
%   a: Vector of hyperparameter values
%   n: Size of the sample to draw
% Outputs:
%   r: Random sample vector of size n drawn from a Dir(a) distribution
function r = drchrnd(a, n)
    p = length(a);
    r = gamrnd(repmat(a, n, 1), 1, n, p);
    r = r ./ repmat(sum(r, 2), 1, p);
end