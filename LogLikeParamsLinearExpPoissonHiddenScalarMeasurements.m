function [logp, p] = ...
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
%
%%

%% Defaults
neuralModel_default = @(p,x)(LinExp(p,x));

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'params')            % Parameters of model
addRequired(Parser,'s')                 % Samples
addRequired(Parser,'o')                 % Production values
addRequired(Parser,'r')                 % M by N (trials by neurons) matrix of neural spike counts
% addParameter(Parser,'neuralModel',neuralModel_default)       % Model of activation functions
addParameter(Parser,'measurementN',2)                        % Number of measurements
addParameter(Parser,'sMinMax',NaN)                           % Min and max of sample interval
addParameter(Parser,'poolOn',false)                          % For controlling parallel pool size; defaults to no pool

parse(Parser,params,s,o,r,varargin{:})

params = Parser.Results.params;
s = Parser.Results.s;
o = Parser.Results.o;
r = Parser.Results.r;
% neuralModel = Parser.Results.neuralModel;
measurementN = Parser.Results.measurementN;
sMinMax = Parser.Results.sMinMax;
poolOn = Parser.Results.poolOn;

% Validate inputs
s = s(:);
o = o(:);
M = length(s);
if size(r,1) == M
    N = size(r,2);
else
    error('Number of enteries in s/o must equal number of rows in r!')
end

if any(isnan(sMinMax))
    smin = min(s);
    smax = max(s);
end

wm = params(1);
wo = params(2);
neuralP = reshape(params(3:end),[N,2]);

%% Integrate out the hidden measurement values and estimate likelihood of model

method_options.dx = 0.01;
functionHandle = @(m,Y)(integrandFunction(m,Y(:,1),Y(:,2),Y(:,3:end),wm,wo,neuralP,smin,smax,poolOn));
p = ndintegrate(functionHandle,...
    repmat([smin-2*wm*smin smax+3*wm*smax],measurementN,1),...
    'method','quad','options',method_options,'ExtraVariables',[s, o, r]);


logp = sum(log(p));
p = exp(logp);              % Typically too small for numerical precision


%% Functions

% Linear exponential function
function lambda = LinExp(params,x)
for triali = 1:size(x,3)
    l = params*[x(:,:,triali) ones(size(x,1),1)]';       % Linear transformation
    lambda(:,triali) = exp( l' );                      % Exponential nonlinearity
end

% Poisson distribution
function p = PoissonDistribution(k,lambda)
p = lambda.^repmat(k',[size(lambda,1),1]) .* exp(-lambda) ./ factorial(repmat(k',[size(lambda,1),1]));

function out =  integrandFunction(m,s,o,r,wm,wo,neuralP,smin,smax,poolOn)
% Set up variables
M = repmat(m,[1, 1, size(s,1)]);
S = repmat( permute(s,[3,2,1]), [size(M,1) 1 1]  );
O = repmat( permute(o,[3,2,1]), [size(M,1) 1 1]  );


% Probability of m given s
if isnan(s)
    probMtakeS = 1;
else
    N = size(m,2);
    probMtakeS = (1/sqrt(2*pi)/wm./S).^N .* exp( -(sum((repmat(S,[1,N,1])-M).^2,2))./(2*wm.^2.*S.^2) );
end

% Probability of o given, m
if isnan(o)
    probOtakeM = 1;
else
    method_opts.dx = 10;%dt;
    method_opts.type = 'quad';
    e = ScalarBayesEstimators(m,wm,smin,smax,'method',method_opts);
    E = repmat(e,[1, 1, size(o,1)]);
    probOtakeM = ( 1/sqrt(2*pi)/wo./E ) .* exp( -(O-E).^2 ./ (2*wo^2*E.^2) );
end

% Probability of r given m
if poolOn
    ptemp = ones(size(probOtakeM));
    ptemp = repmat(ptemp,[1,1,1,size(r,2)]);
    parfor neuroni = 1:size(r,2)
        %     LAMBDA = LinExp(neuralP(neuroni,:),M);
        LAMBDA = LinExp(neuralP(neuroni,:),M(:,1,:));
        ptemp(:,:,:,neuroni) = permute(...
            PoissonDistribution(r(:,neuroni),LAMBDA),...
            [1 3 2]);
    end
    probRtakeM = prod(ptemp,4);
else
    probRtakeM = ones(size(probOtakeM));
    for neuroni = 1:size(r,2)
        %     LAMBDA = LinExp(neuralP(neuroni,:),M);
        LAMBDA = LinExp(neuralP(neuroni,:),M(:,1,:));
        ptemp = PoissonDistribution(r(:,neuroni),LAMBDA);
        probRtakeM = probRtakeM .* permute(ptemp,[1 3 2]);
    end
end

out = permute( probMtakeS .* probOtakeM .* probRtakeM, [1,3,2]);