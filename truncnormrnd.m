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
    if (nargin < 2) || isempty(mu)      % Default mean
      mu = 0;
    end
    if (nargin < 3) || isempty(sig)     % Default standard deviation
      sig = 1;
    end
    if (nargin < 4) || isempty(xlo)     % Default low truncation point
      xlo = -inf;
      plo = 0;
    else
      plo = normcdf((xlo - mu) / sig); % Standard normal percentile for low truncation point
    end
    if (nargin < 5) || isempty(xhi)     % Default high truncation point
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