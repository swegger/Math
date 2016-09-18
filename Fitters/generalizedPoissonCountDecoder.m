function [e, logLikelihood, feval, exitflg, output] =...
    generalizedPoissonCountDecoder(r,theta,x_init,varargin)
%% tsSpikeCountDecoder
%
%   e = generalizedPoissonCountDecoder(r,theta,x_init)
%   Estimates, trial by trial, some unknonw quantity represent by a set of
%   firing rates, r, with linear tuning parameters theta. Starts
%   optimiztion from initial guess, x_init.
%
%   e = generalizedPoissonCountDecoder(...,'activationFunction',fun)
%   Allows user to supply own activation function that relates the
%   dependent variable to the mean count output.
%
%%

%% Defaults
OPTIONS = optimset('Display','iter');
defaultFunction = @(x,theta)(theta(1)*x + theta(2));

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'r');                                        % m x n Matrix of firing rates with m trials and n units
addRequired(Parser,'theta');                                    % 2 x n matrix of tuning parameters
addRequired(Parser,'x_init',NaN);                               % Initial estimate
addParameter(Parser,'options',OPTIONS)                          % Options for fminsearch
addParameter(Parser,'activationFunction',defaultFunction)       % Relationship between dependent variable and mean count

parse(Parser,r,theta,x_init,varargin{:})

r = Parser.Results.r;
theta = Parser.Results.theta;
x_init = Parser.Results.ts_init;
options = Parser.Results.options;
activationFunction = Parser.Results.activationFunction;

%% Set up gradient descent algorithm
for i = 1:size(r,1)
    f = @(ts)(-ML(ts,r(i,:),theta,activationFunction));
    [e(i), feval(i), exitflg(i), output(i)] = fminsearch(f, x_init(i), options);
    logLikelihood(i) = -feval(i);
end

%% Functions
function out = ML(x,r,theta,activationFunction)
llike = r .*...
    log( activationFunction(x,theta) ) - ( activationFunction(x,theta) );

out = sum(llike,2);
