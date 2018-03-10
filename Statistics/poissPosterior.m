function [p, lambdaML, CDF] = poissPosterior(k,lambda)
%% poissPosterior
%
%   p = poissPosterior
%
%   Finds the probability of a given lambda given a set of measured counts,
%   k.
%
%%

%% Defaults
options_default.dx = 0.001;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'k')
addRequired(Parser,'lambda')
addParameter(Parser,'options',options_default)

parse(Parser,k,lambda)

k = Parser.Results.k;
lambda = Parser.Results.lambda;
options = Parser.Results.options;

%% Calc posterior
Likelihood = numerator(k,lambda);
f = @(lambda)(numerator(k,lambda));
Z = ndintegrate(f,[0 mean(k)*10],'method','quad','options',options);
p = 1/Z * Likelihood * (lambda(2) - lambda(1));

%% Cumulative distribution function
% g = @(lambda,k)(
% CDF = ndintegrate(g,[0 mean(k)*10]
CDF = [];

%% Maximum likelihood lambda
lambdaML = mean(k);

%% Functions

%% Numerator
function likelihood = numerator(k,lambda)
%%
    LAMBDA = repmat(lambda(:),[1, numel(k)]);
    K = repmat(k(:)',[numel(lambda), 1]);
    L = LAMBDA.^(K).*exp(-LAMBDA)./factorial(K);
    likelihood = prod(L,2);

% 