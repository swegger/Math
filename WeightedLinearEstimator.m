function [e, likelihood] = WeightedLinearEstimator(m,sigM,sigP,muP,varargin)
%% WeightedLinearEstimator
%
%   [e, likelihood] = WeigthedLinearEstimator(m,sigM,sigP,muP)
%   Finds the estimate, e, of an input given the measurement(s), m under a
%   Bayesian model that assumes Gaussian measurement noise model with
%   variance sigM and a Gaussian prior with variance sigP and mean muP.
%
%
%%

%% Defaults
method_opts.type = 'integral';
estimator_opts.type = 'BLS';
prior_opts.type = 'uniform';

%% Parse inputs
inputs = inputParser;
addRequired(inputs,'m');
addRequired(inputs,'sigM');
addRequired(inputs,'sigP');
addRequired(inputs,'muP');

parse(inputs,m,sigM,sigP,muP,varargin{:})

m = inputs.Results.m;
sigM = inputs.Results.sigM;
sigP = inputs.Results.sigP;
muP = inputs.Results.muP;


%% Compute the estimate

% Number of measurements
N = size(m,2);
if length(sigM) == 1
    sigM = sigM*ones(1,N);
elseif length(sigM) ~= N
    error('Number of measurement noise parameters must be one or equal to number of measurements!')
end

% Determine weights
wP = (1/sigP) / ( sum(1./sigM) + 1/sigP );
for i = 1:N
    w(1,i) = (1/sigM(i)) / ( sum(1./sigM) + 1/sigP );
end

% Compute estimate
e = sum( repmat(w,[size(m,1),1]) .* m ,2) + wP*muP;

