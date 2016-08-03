function theta = FitLDS_PoissonObs(s,N,theta0)
%% FitLDS_PoissonObs
%
%
%%

%% Defaults

%% Parse input
Parser = inputParser;

addRequired(Parser,'s')
addRequired(Parser,'N')
addRequired(Parser,'theta0')

parse(Parser,s,N)

s = Parser.Results.s;
N = Parser.Results.N;
theta0 = Parser.Results.theta0;

%% Fit the model to the data using EM

% Compute E-step
E = expectationStep(y,

% Compute M-step


%% Functions
function E = expectationStep(s,N,theta);

