function logp = ...
    LogLikeParamsLinearExpPoissonHiddenScalarMeasurements(...
        params,s,o,r,varargin)
%% LogLikeParamsLinearExpPoissonHiddenScalarMeasurements
%
%   logp = LogLikeParamsLinearExpPoissonHiddenScalarMeasurements(...
%        params,s,o,r)
%
%       Determines log-likelihood of parameters of
%       linear-exponential-poisson neural model with dependence on hidden
%       measurements with scalar variability given sample inputs (s),
%       behavioral outputs, o, and neural responses, r.
%
%%

%% Defaults


%% Parse inputs
Parser = inputParser;

addRequired(Parser,'params')            % Parameters of model
addRequired(Parser,'s')                 % Samples
addRequired(Parser,'o')                 % Production values
addRequired(Parser,'r')                 % M by N (trials by neurons) matrix of neural spike counts

parse(Parser,params,s,o,r,varargin{:})

params = Parser.Results.params;
s = Parser.Results.s;
o = Parser.Results.o;
r = Parser.Results.r;

% Validate inputs
s = s(:);
o = o(:);
M = length(s);
if size(r,1) == M
    N = size(r,2);
else
    error(['Number of enteries in s/o must equal number of rows in r!')
end

%% Do stuff

%% Functions

% Scalar noise model
function p = ScalarNoise(w,x,y)
p = 1./sqrt(2*pi*w^2*x^2) .* exp( -(y-x).^2./(2*w^2*x^2) );

% Sca