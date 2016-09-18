function [u, FVAL, EXITFLAG, OUTPUT] = minimizeRayleighQuotient(u0,A,B,varargin)
%% minimizeRayleighQuotient
%
%   u = minimizeRayleighQuotient(A,B)
%
%   Finds the vector, u, that minimizes the Rayleigh quotient of the form:
%      ( u'*A*u ) / ( u'*B*u)
%
%%

%% Defaults
OPTIONS_default = optimset('Display','iter');

%% Parse inputs

Parser = inputParser;
addRequired(Parser,'u0')
addRequired(Parser,'A')
addRequired(Parser,'B')
addParameter(Parser,'OPTIONS',OPTIONS_default)

parse(Parser,u0,A,B,varargin{:})

u0 = Parser.Results.u0;
A = Parser.Results.A;
B = Parser.Results.B;
OPTIONS = Parser.Results.OPTIONS;

%% Contruct minimization function
f = @(u)( RayleighQuotient(u,A,B) );

%% Perform minimization
[u, FVAL, EXITFLAG, OUTPUT] = fminsearch(f,u0,OPTIONS);


%% Functions
function quotient = RayleighQuotient(u,A,B)
quotient = (u'*A*u) / (u'*B*u);