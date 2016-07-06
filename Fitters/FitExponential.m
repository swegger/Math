function [theta, feval, exitflg, output] = FitExponential(x,y,theta_init,varargin)
%% FitGaussianMAP
%
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
f = @(p)(objective(p,x,y));
[theta, feval, exitflg, output] = fmincon(f, theta_init, [],[],[],[],[0 -100],[0.01 100],[],options);

%% Functions
% Log posterior of the parameters, given the data
function out = objective(p,x,y)
logyhat = p(1)*(x - p(2));
error = (log(y) - logyhat).^2;
out = sum(error);