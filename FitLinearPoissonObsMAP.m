function [theta, logPosterior, exitflg, output] = FitLinearPoissonObsMAP(x,y,theta_init,varargin)
%% FitGaussianMAP
%
%   theta = FitLinearPoissonObsMAP(x,y,theta_init)
%
%   Fits a linear model to data with Poisson observation noise. Model assumes
%   the following form:
%       p(y|x) = (lambda^y*exp(-lambda))/(y!)
%       lambda = theta(1)*x + theta(2);
%   By default, a prior over the parameters, theta, is not considered; thus
%   the alorithm finds the maximum likelihood estimate of the parameter
%   values.
%
%   [theta, logLikelihood] = FitLinearPoissonObsMAP(x,y,theta_init)
%   Returns the value of the log-likelihood at the best fitting parameters.
%
%   [theta, logPosterior] =
%   FitLinearPoissonObsMAP(x,y,theta_init,'mu',mu,'sig',sig)
%   Fits a linear model to the data, assuming the prior distribution of
%   parameter values is a multivariate Gaussian centered on mu and with
%   covariance sig.
%
%   
% written by swe 20160412
%%

%% Defaults
OPTIONS = optimset('Display','iter');

%% Parse inputs
Parser = inputParser;
addRequired(Parser,'x');                        % Dependent variable
addRequired(Parser,'y');                        % Observations
addRequired(Parser,'theta_init');               % Initial values
addParameter(Parser,'mu',NaN);                  % Mean of prior
addParameter(Parser,'sig',NaN);                 % Variance of prior
addParameter(Parser,'options',OPTIONS)          % Options for fminsearch

parse(Parser,x,y,theta_init,varargin{:})

x = Parser.Results.x;
y = Parser.Results.y;
theta_init = Parser.Results.theta_init;
mu = Parser.Results.mu;
sig = Parser.Results.sig;
options = Parser.Results.options;

if isnan(mu)
    mu = zeros(size(theta_init));
    mu = mu(:);
end
if isnan(sig)
    sig = eye(length(theta_init));
    sig(sig~=0) = Inf;
end

%% Fit the data
f = @(p)(-MAP(p,x,y,mu,sig));
[theta, feval, exitflg, output] = fminsearch(f, theta_init, options);
logPosterior = -feval;

%% Functions
% Log posterior of the parameters, given the data
function out = MAP(p,x,y,mu,sig)
p = p(:);
lambda = @(x)(p(1)*x + p(2));
loglikelihood = y.*log(lambda(x)) - lambda(x);
warning('off','all')
logprior = -1/2 * (p-mu)'*inv(sig)*(p-mu);
warning('on','all')
logposterior = loglikelihood + logprior;
out = sum(logposterior);