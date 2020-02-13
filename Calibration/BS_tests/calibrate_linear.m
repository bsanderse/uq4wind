function [x, w, err, n, model_evals] = calibrate_linear(u, varargin)
% [x, w, err, n] = CALIBRATE(u, ...) - Construct calibration case and pass
% it to analyticnd.
%
% This function constructs a calibration case of the following basic model:
%   z = u(x0) + eps,
% where eps is normally distributed using known standard deviation.
%
% u should be a function that accepts a Nxd matrix (d can be changed, see
% below) and returns a Nx1 vector that represents the function outputs at
% all locations. This function assumes u can be evaluated fast, as it does
% not cache intermediate results.
%
% All options passed to this function are passed to analyticnd, except for
% the options gobbled by this function (all are optional):
% nz    - Number of data points that should be used
%         Default: 20
% sigma - Standard deviation of eps
%         Default: 1/10
% x0    - Point where artificial data is produced (i.e. z = u(x0))
%         Default: 1/2
% d     - Dimension of the problem
%         Default: 5
% prior - If true, does not uses the posterior, but the prior for the
%         calibration. Can be used for comparison.
%         Default: false
% 
% truemean/I are determined by this function and cannot be overwritten.
%
% Output of the function is exactly the output obtained by running
% analyticnd. See the documentation of that function for more information
% about that.
%
% See also ANALYTICND

%% Parse optional function arguments
% nz = 20;
sigma = 1/10;
d = 5;
% x0 = [];
S = 1e5; % Used for estimation of truemean
use_prior = false;
varargin_analyticnd = {};
for i=1:2:length(varargin)
    switch varargin{i}
        case 'nz'
            nz = varargin{i+1};
        case 'sigma'
            sigma = varargin{i+1};
        case 'x0'
            x0 = varargin{i+1};
        case 'd'
            d = varargin{i+1};
        case 'prior'
            use_prior = varargin{i+1};
        case 'S'
            S = varargin{i+1}*100;
            varargin_analyticnd{end+1} = varargin{i};
            varargin_analyticnd{end+1} = varargin{i+1};
        case 'A'
            A = varargin{i+1};
        case {'truemean', 'I'} % Determined in this function.
            warning('Option %s ignored.', varargin{i});
        case 'z'
           z = varargin{i+1}; 
        otherwise
            varargin_analyticnd{end+1} = varargin{i};
            varargin_analyticnd{end+1} = varargin{i+1};
    end
end

% if isempty(x0)
%     x0 = ones(1, d)/2;
% end
% 
% if length(x0) ~= d
%     warning('Size of x0 does not match d. Resetting to default.');
%     x0 = ones(1, d)/2;
% end

%% Construct statistical model, so data, prior, likelihood, posterior
I = [0, 1];
% z = u(x0,A) + normrnd(0, sigma, 1, nz);
prior = @(x) ones(size(x, 1), 1);
% 1/ sqrt( (2*pi)^d * det(Sigma) ) * exp ( -(1/2) * (u-z)' * Sigma^{-1} *
% (u-z) )
% for the simple case Sigma = sigma^2*I, we have
% det(Sigma) = (sigma^2)^d
% Sigma^{-1} = 1/sigma^2 * I
% this gives the multiplicative factor 1/ sqrt( (2*pi)^d * (sigma^2)^d ); 
% however, this is just a scaling that is not important for the analysis
likelihood = @(x) exp(- sum((u(x,A)-z).^2, 2) ./ (2*sigma^2) ); 
posterior = @(x) likelihood(x).*prior(x);

%% Construct (approximation) of true mean using MC
% X = rand(S, d); % this assumes uniform distributed RV
X = rand(1e5,d);
post     = posterior(X);
evidence = mean(post); % note, this involves S model evaluations
truemean = mean(X.*post) / evidence;
clear X

%% Finally, overwrite sampler if user wants that
if use_prior
    varargin_analyticnd{end+1} = 'sampler';
    varargin_analyticnd{end+1} = prior;
    
    % Scale posterior with evidence
    posterior = @(x) likelihood(x).*prior(x) / evidence;
end

%% Run analyticnd with the obtained posterior and truemean estimate
[x, w, err, n, model_evals] = analyticnd(d, posterior, ...
    'truemean', truemean, ...
    'I', I, ...
    varargin_analyticnd{:});

end
