function [e, logLikelihood, exitflg, output] = LinearExponentialSpikeCountDecoder(r,theta,varargin)
%% tsSpikeCountDecoder
%
%   e = LinearExponentialSpikeCountDecoder(r,theta)
%   Estimates a sample, s, trial by trial, from a set of firing rates, r,
%   with linear tuning parameters theta.
%
%%

%% Defaults
OPTIONS = optimset('Display','iter');

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'r');                    % m x n Matrix of firing rates with m trials and n units
addRequired(Parser,'theta');                % 2 x n matrix of tuning parameters
addParameter(Parser,'ts_init',NaN);         % Initial conditions for gradient descent estimator
addParameter(Parser,'options',OPTIONS)          % Options for fminsearch
addParameter(Parser,'estimator','fminsearch')   % Estimation algorithm

parse(Parser,r,theta,varargin{:})

r = Parser.Results.r;
theta = Parser.Results.theta;
ts_init = Parser.Results.ts_init;
options = Parser.Results.options;
estimator = Parser.Results.estimator;

if isnan(ts_init)
    ts_init = 800*ones(size(r,1),1);
end

%% ML estimate of s
switch estimator
    case 'SumLogApprox'
        error('SumLogApprox algorithm is not yet functional')
        thetas = repmat( permute(theta,[3,2,1]) , [size(r,1),1,1]);
        rs = repmat(r,[1,1,size(theta,1)]);
        numerator = r.*thetas(:,:,1) - log(thetas(:,:,1)) - thetas(:,:,2);
        e = sum(numerator,2)./sum(thetas(:,:,1),2);
        logLikelihood = ML(e,r,theta);
        exitflg = NaN;
        output = NaN;
        
    case 'fminsearch'
        % Set up gradient descent algorithm
        for i = 1:size(r,1)
            f = @(ts)(-ML(ts,r(i,:),theta));
            %[te(i), feval(i), exitflg(i), output(i)] = fminsearch(f, ts_init(i), options);
            [e(i), feval(i), exitflg(i), output(i)] = fmincon(f, ts_init(i), [], [], [], [], [0],[Inf], [], options);
            logLikelihood(i) = -feval(i);
        end
        
    otherwise
        error(['Estimator type ' estimator ' not recognized!'])
end


%% Functions
function out = ML(ts,r,theta)
rhat = exp( repmat(theta(1,:),[size(r,1),1]).*repmat(ts,[1,size(theta,2)]) + repmat(theta(2,:),[size(r,1),1]) );
rhat(rhat < 0) = 0;
llike = r .*...
    log( rhat )...
    - ( rhat );
%llike( repmat(theta(1,:),[size(r,1),1]).*ts + repmat(theta(2,:),[size(r,1),1]) < 0 ) = -10000;

out = nansum(llike,2);
