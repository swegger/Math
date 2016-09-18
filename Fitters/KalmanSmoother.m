function [mu, V] = kalmanSmoother(x,A,C,Gamma,Sigma,mu0,V0,varargin)
%% kalmanSmoother
%
%   [mu, V] = kalmanSmoother(x,A,C,Gamma,Sigma,mu0,V0)
%
%   Computes the optimal estimate of the hidden state at each time point,
%   z(k), using the full set of observations x.
%
%%

%% Defaults

%% Parse input
Parser = inputParser;

addRequired(Parser,'x')
addRequired(Parser,'A')
addRequired(Parser,'C')
addRequired(Parser,'Gamma')
addRequired(Parser,'Sigma')
addRequired(Parser,'mu0')
addRequired(Parser,'V0')
addParameter(Parser,'K',Inf)

parse(Parser,x,A,C,Gamma,Sigma,mu0,V0,varargin{:})

x = Parser.Results.X;
A = Parser.Results.A;
C = Parser.Results.C;
Gamma = Parser.Results.Gamma;
Sigma = Parser.Results.Sigma;
mu0 = Parser.Results.mu0;
V0 = Parser.Results.V0;
K = Parser.Results.K;

if isinf(K)
    K = size(x,2);
end

hiddenDims = size(A,1);

%% Perform smoothing
