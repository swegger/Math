function [logp, p] = ScalarMeasurementProbability(...
    m,wm,smin,smax,varargin)
%% ScalarMeasurementProbabilityGivenInputOutput
%
%   p = ScalarMeasurementProbabilityGivenInputOutput(s,o,m,wm,wp,smin,smax)
%       Determines the probability of a vector of n measurements given the
%       input, s, to a system and its output, o, under the BLS model for
%       optimal estimation with a uniform prior distribution between smin
%       and smax.
%
%%

%% Defaults
EstimatorOptions_default.method_opts.type = 'quad';
EstimatorOptions_default.method_opts.dx = 10;
EstimatorOptions_default.estimator_opts.type = 'BLS';
EstimatorOptions_default.prior.type = 'uniform';
method_options_default.dx = 0.01;
Normalization_default = false;

%% Parse input
Parser = inputParser;

addRequired(Parser,'m')
addRequired(Parser,'wm')
addRequired(Parser,'smin')
addRequired(Parser,'smax')
addParameter(Parser,'method_options',method_options_default)
addParameter(Parser,'Normalization',Normalization_default)

parse(Parser,m,wm,smin,smax,varargin{:})

m = Parser.Results.m;
wm = Parser.Results.wm;
smin = Parser.Results.smin;
smax = Parser.Results.smax;
method_options = Parser.Results.method_options;
Normalization = Parser.Results.Normalization;


%% Probability of m
functionHandle = @(s,Y)(integrandFunction(s,Y,wm));
p = ndintegrate(functionHandle,[smin smax],'method','quad',...
    'options',method_options,'ExtraVariables',m);

logp = log(p);

%% Normalization factor
% switch Normalization
%     case {'Yes','yes','Y','y',true,1}
%         warning('Unclear if normalization is working properly')
%         method_options.dx = 0.01;
%         functionHandle = @(m,Y)(integrandFunction(m,Y(:,1),Y(:,2),wm,wo,smin,smax));
%         Z = ndintegrate(functionHandle,repmat([smin-2*wm*smin smax+3*wm*smax],N,1),'method','quad','options',method_options,'ExtraVariables',[s, o]);
%     case {'No','no','N','n',false,0}
%         Z = ones(1,numel(s));
%     otherwise
%         error('Normalization option not recognized!')
% end

%% Functions
function out =  integrandFunction(s,m,wm)
% Set up variables
N = size(m,2);
M = repmat(m,[1, 1, size(s,1)]);
S = repmat( permute(s,[3,2,1]), [size(M,1) 1 1]  );

% Probability of m given s
probMtakeS = (1/sqrt(2*pi)/wm./S).^N .* exp( -(sum((repmat(S,[1,N,1])-M).^2,2))./(2*wm.^2.*S.^2) );

% out = permute( probS .* probMtakeS .* probOtakeM , [1,3,2]);  %p(s) shouldn't appear here...
out = permute( probMtakeS , [3,1,2]);