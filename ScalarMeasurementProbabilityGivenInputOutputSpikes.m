function [logp, p] = ScalarMeasurementProbabilityGivenInputOutputSpikes(...
    s,o,r,neuralModel,m,wm,wo,smin,smax,varargin)
%% ScalarMeasurementProbabilityGivenInputOutput
%
%   p = InferMeasurementsWithSpikes(s,o,r,m,wm,wp,smin,smax)
%       Determines the probability of a vector of n measurements given the
%       input, s, to a system and its output, o, under the BLS model for
%       optimal estimation with a uniform prior distribution between smin
%       and smax and a spike enocoding model of the form:
%
%           p(r_i|m) = \lambda_i^{r_i}exp(-\lambda_i) / r_i!
%
%       with \lambda_i being the deterministic encoding model between m_1, 
%       m_2 and the response of neuron i, r_i. The form of the encoding
%       model and its parameters should be stored in the structure
%       neuralModel with the fields:
%           neuralModel(i).f - user defined function of m
%                              neuralModel(i).params
%                            - e.g. f(m,neuralModel(i).params)
%           neuralModel(i).params - parameters of model fit to data
%           
%
%%

%% Defaults
EstimatorOptions_default.method_opts.type = 'quad';
EstimatorOptions_default.method_opts.dx = 10;
EstimatorOptions_default.estimator_opts.type = 'BLS';
EstimatorOptions_default.prior.type = 'uniform';
Normalization_default = false;

%% Parse input
Parser = inputParser;

addRequired(Parser,'s')
addRequired(Parser,'o')
addRequired(Parser,'r')
addRequired(Parser,'neuralModel')
addRequired(Parser,'m')
addRequired(Parser,'wm')
addRequired(Parser,'wo')
addRequired(Parser,'smin')
addRequired(Parser,'smax')
addParameter(Parser,'EstimatorOptions',EstimatorOptions_default)
addParameter(Parser,'Normalization',Normalization_default)

parse(Parser,s,o,r,neuralModel,m,wm,wo,smin,smax,varargin{:})

s = Parser.Results.s;
o = Parser.Results.o;
r = Parser.Results.r;
neuralModel = Parser.Results.neuralModel;
m = Parser.Results.m;
wm = Parser.Results.wm;
wo = Parser.Results.wo;
smin = Parser.Results.smin;
smax = Parser.Results.smax;
EstimatorOptions = Parser.Results.EstimatorOptions;
Normalization = Parser.Results.Normalization;


% Ensure s and o are the proper dimensions
if any(isnan(s)) && numel(s) ~= numel(o)
    s = nan(size(o));
elseif any(isnan(o)) && numel(s) ~= numel(o)
    o = nan(size(s));
elseif numel(s) ~= numel(o)
    error('Input, s, and output, o, must have equal dimensions!')
end
s = s(:);
o = o(:);
N = size(m,2);


%% Expand s, m and o to proper dimenions
M = repmat(m,[1, 1, size(s,1)]);
S = repmat( permute(s,[3,2,1]), [size(M,1) 1 1]  );
O = repmat( permute(o,[3,2,1]), [size(M,1) 1 1]  );

%% Probability of s
if isnan(s)
    probS = 1;
else
    probS = 1 ./ (smax - smin);     % Assumes uniform distribution and is therefore independent of the actual s
end

%% Probability of m given s
if isnan(s)
    probMtakeS = 1;
else
    probMtakeS = (1/sqrt(2*pi)/wm./S).^N .* exp( -(sum((repmat(S,[1,N,1])-M).^2,2))./(2*wm.^2.*S.^2) );
end

%% Probability of o given m
if isnan(o)
    probOtakeM = 1;
else
    e = ScalarBayesEstimators(m,wm,smin,smax,'method',EstimatorOptions.method_opts,...
        'estimator',EstimatorOptions.estimator_opts,'prior',EstimatorOptions.prior);
    E = repmat(e,[1, 1, size(o,1)]);
    probOtakeM = ( 1/sqrt(2*pi)/wo./E ) .* exp( -(O-E).^2 ./ (2*wo^2*E.^2) );
end

%% Probability of r given m
logprobRtakeM = zeros(size(probOtakeM));
for neuroni = 1:size(r,2)
    lambda = repmat(neuralModel(neuroni).f(M,neuralModel(neuroni).params),...
        [1, 1, size(s,1)]);
    probabilisticResponse = log(...
        lambda.^repmat(permute(r(:,neuroni),[3 2 1]),[size(lambda,1), 1, 1])...
        .* exp( -lambda ) ./ ...
        factorial(repmat(permute(r(:,neuroni),[3 2 1]),[size(lambda,1), 1, 1]))...
        );
                          
    logprobRtakeM = logprobRtakeM + probabilisticResponse;
end

%% Normalization factor
switch Normalization
    case {'Yes','yes','Y','y',true,1}
        warning('Unclear if normalization is working properly')
        method_options.dx = 0.01;
        functionHandle = @(m,Y)(integrandFunction(m,Y(:,1),Y(:,2),wm,wo,smin,smax));
        Z = ndintegrate(functionHandle,repmat([smin-2*wm*smin smax+3*wm*smax],N,1),'method','quad','options',method_options,'ExtraVariables',[s, o]);
    case {'No','no','N','n',false,0}
        Z = ones(1,numel(s));
    otherwise
        error('Normalization option not recognized!')
end

%% Combine all the evidence
% logp = permute(...
%     log(probS) + log(probMtakeS) + log(probOtakeM) + logprobRtakeM,...
%     [1,3,2]) - ...
%     repmat(Z,size(m,1),1);        % p(s) shouldn't appear here...

logp = permute(...
    log(probMtakeS) + log(probOtakeM) + logprobRtakeM,...
    [1,3,2]) - ...
    repmat(Z,size(m,1),1);

p = exp(logp);


%% Functions
function out =  integrandFunction(m,s,o,wm,wo,smin,smax)
% Set up variables
M = repmat(m,[1, 1, size(s,1)]);
S = repmat( permute(s,[3,2,1]), [size(M,1) 1 1]  );
O = repmat( permute(o,[3,2,1]), [size(M,1) 1 1]  );

% Probability of s
if isnan(s)
    probS = 1;
else
    probS = 1 ./ (smax - smin);     % Assumes uniform distribution and is therefore independent of the actual s
end;

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

% out = permute( probS .* probMtakeS .* probOtakeM , [1,3,2]);  %p(s) shouldn't appear here...
out = permute( probMtakeS .* probOtakeM , [1,3,2]);