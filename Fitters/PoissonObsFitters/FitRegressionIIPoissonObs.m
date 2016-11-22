function [theta, logPosterior, exitflg, output] = FitRegressionIIPoissonObs(x,y,theta_init,varargin)
%% FitGaussianMAP
%
%   theta = FitRegressionIIPoissonObs(x,y,theta_init)
%
%   Fits a linear model to data with Poisson observation noise. Model assumes
%   the following form:
%       p(y|x) = (lambda^y*exp(-lambda))/(y!)
%       lambda = theta(1)*x + theta(2);
%   By default, a prior over the parameters, theta, is not considered; thus
%   the alorithm finds the maximum likelihood estimate of the parameter
%   values.
%
%   [theta, logLikelihood] = FitRegressionIIPoissonObs(x,y,theta_init)
%   Returns the value of the log-likelihood at the best fitting parameters.%
%   
% written by swe 20160415
%%

%% Defaults
OPTIONS = optimset('Display','iter');

%% Parse inputs
Parser = inputParser;
addRequired(Parser,'x');                        % Dependent variable
addRequired(Parser,'y');                        % Observations
addRequired(Parser,'theta_init');               % Initial values
addParameter(Parser,'options',OPTIONS)          % Options for fminsearch

parse(Parser,x,y,theta_init,varargin{:})

x = Parser.Results.x;
y = Parser.Results.y;
theta_init = Parser.Results.theta_init;
options = Parser.Results.options;

%% Fit the data
f = @(p)(-ML(p,x,y));
[theta, feval, exitflg, output] = fminsearch(f, theta_init, options);
logPosterior = -feval;

%% Functions
% Log posterior of the parameters, given the data
function out = ML(p,x,y)
p = p(:);
logLikelihood = x*log(p(1)) - p(1) + y.*log(x*p(2)+p(3)) - x*p(2) - p(3);
out = sum(logLikelihood);