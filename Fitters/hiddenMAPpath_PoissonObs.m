function x = hiddenMAPpath_PoissonObs(s,theta,t)
%% hiddenMAPpath_PoissonObs
%
%
%
%%

%% Defaults

%% Parse input
Parser = inputParser;

addRequired(Parser,'s')
addRequired(Parser,'theta')
addRequired(parser,'t')

parse(Parser,s,theta,t)

s = Parser.Results.s;
theta = Parser.Results.theta;
t = Parser.Results.t;

