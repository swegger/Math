function [theta, logPosterior, exitflg, output] = FitConstantPoissonObsMAP(x,y,varargin)
%% FitGaussianMAP
%
%   theta = FitConstantPoissonObsMAP(x,y,theta_init)
%
%   Fits a linear model to data with Poisson observation noise. Model assumes
%   the following form:
%       p(y|x) = (lambda^y*exp(-lambda))/(y!)
%       lambda = theta;
%   By default, a prior over the parameters, theta, is not considered; thus
%   the alorithm finds the maximum likelihood estimate of the parameter
%   values using the analytical solution.
%
%   theta = FitConstantPoissonObsMAP(...,'algorithm',algorithm)
%   Allows user to choose the algorithm to solve for the MAP solution.
%       'analytical' - use the analytical solution (default)
%       'fminsearch' - use fminsearch (not recomended)
%
%   theta = FitConstantPoissonObsMAP(...,'prior',prior)
%   Allows user to specify the prior type.
%       'Gamma' - use a Gamma distribution (default)
%       'Gaussian' - use a Gaussian; note, analytical solution is not known
%
%   [theta, logLikelihood] = FitConstantPoissonObsMAP(x,y,theta_init)
%   Returns the value of the log-likelihood at the best fitting parameters.
%
%   [theta, logPosterior] =
%   FitConstantPoissonObsMAP(x,y,theta_init,'a',a,'b',b)
%   Fits a linear model to the data, assuming the prior distribution of
%   parameter values. If the prior is a Gamma distribution (default), a and
%   b parameterize the distirbution as:
%           p(theta) = (1/gamma(a)) * b^a * theta^(a-1) * exp(-b*theta)
%   If the prior is a multivariate Gaussian, a and b parameterize the
%   distribution as:
%           p(theta) = (1/sqrt(2*pi*b^2)) * exp( -(theta-a)^2 / (2*b^2) )
%   Generally, a Gaussian prior is not recommended because it allows
%   negative values of theta, which does not make sense in the context of a
%   Poisson measurement system.
%
%
%   Written by swe 20160412
%%

%% Defaults
OPTIONS = optimset('Display','iter');

%% Parse inputs
Parser = inputParser;
addRequired(Parser,'x');                        % Dependent variable
addRequired(Parser,'y');                        % Observations
addParameter(Parser,'theta_init',NaN);          % Initial values
addParameter(Parser,'a',NaN);                   % Mean of prior
addParameter(Parser,'b',NaN);                   % Variance of prior
addParameter(Parser,'options',OPTIONS)          % Options for fminsearch
addParameter(Parser,'algorithm','analytical')   % Algorithm to use
addParameter(Parser,'prior','Gamma')            % Prior to implement

parse(Parser,x,y,varargin{:})

x = Parser.Results.x;
y = Parser.Results.y;
theta_init = Parser.Results.theta_init;
a = Parser.Results.a;
b = Parser.Results.b;
options = Parser.Results.options;
algorithm = Parser.Results.algorithm;
prior = Parser.Results.prior;

if strcmp(algorithm,'fminsearch') && isnan(theta_init)
    error('User must supply the initial conditions if one wishes to use fminsearch to find the MAP solution!')
end

if isnan(a)
    switch prior
        case 'Gamma'
            a = 1;
        case 'Gaussian'
            a = 0;
        otherwise
            error([prior ' is not a recognized prior type!']);
    end
end
if isnan(b)
    switch prior
        case 'Gamma'
            b = 0;
        case 'Gaussian'
            b = Inf;
        otherwise
            error([prior ' is not a recognized prior type!']);
    end
end

if strcmp(algorithm,'analytical') && strcmp(prior,'Gaussian')
    warning('Analytical solution for a Gaussian prior is not known; using fminsearch')
    algorithm = 'fminsearch';
end

%% Fit the data
switch algorithm
    case 'analytical'
        switch prior
            case 'Gamma'
                theta = ( sum(y(:)) + (a-1) ) / ( numel(y)+b );
                logPosterior = MAP(theta,x,y,a,b,prior);
                
            case 'Gaussian'
                error('Analytical solution for a Gaussian prior is not known! Choose fminsearch as your algorithm')
                
            otherwise
                error([prior ' is not a recognized prior type!']);
        end
        
    case 'fminsearch'
        f = @(p)(-MAP(p,x,y,a,b,prior));
        [theta, feval, exitflg, output] = fminsearch(f, theta_init, options);
        logPosterior = -feval;
end

%% Functions
% Log posterior of the parameters, given the data
function out = MAP(p,x,y,a,b,prior)
p = p(:);
lambda = @(x)(p(1));
loglikelihood = y.*log(lambda(x)) - lambda(x);
warning('off','all')
switch prior
    case 'Gamma'
        logprior = (a-1)*log(p) - b*p;                % Gamma prior
    case 'Gaussian'
        logprior = -1/2 * (p-a).^2/(b^2);      % Gaussian prior
end
warning('on','all')
logposterior = loglikelihood + logprior;
out = sum(logposterior);