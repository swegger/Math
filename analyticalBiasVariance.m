function [bias, variance] = analyticalBiasVariance(ts,wm,wp,N,varargin)
%% analyticalBiasVariance
%
%   [bias variance] = analyticalBiasVariance(ts,wm,wp,N)
%       Finds mean-squared bias and mean variance of an idea observer using
%       the BLS strategy to estimate timing with N measurements.
%
%   [bias variance] = analytical BiasVariance(...,'method',method)
%       Specifies the method of integration to be used.
%           method = 'quad' - Simpson's quadrature (default)
%           method = 'integral' - TODO;
%           method = 'trapz' - TODO;
%
%   [bias variance] =
%   analyticalBiasVariance(...,'estimator',estimator_opts)
%       Specifies the type of estimator and allow user-defined options
%           estimator_opts.type = 'BLS' - BLS estimator (default)
%           TODO
%
%   [bias variance] = analyticalBiasVariance(...,'options',options)
%       Options to be passed to ta_expectation2
%
%%

% Defaults
method_opts_defaults.dx = 1;
estimator_defaults.type = 'BLS';

% Parse inputs
p = inputParser;
addRequired(p,'ts');
addRequired(p,'wm');
addRequired(p,'wp');
addRequired(p,'N');
addParamValue(p,'method','quad');
addParamValue(p,'method_options',method_opts_defualts);
addParamValue(p,'estimator',estimator_defaults);
addParamValue(p,'options');


parse(p,ts,wm,wp,varargin{:})

ts = p.Results.ts;
wm = p.Results.wm;
wp = p.Results.wp;
N = p.Results.N;
method = p.Results.method;
method_options = p.Results.method_options;
estimator = p.Results.estimator;
options = p.Results.options;

% Constants
MinMax = repmat([min(ts)-3*wm*min(ts) max(ts)+3*wm*max(ts)],N,1);
dt = method_options.dx*(MinMax(1,2)-MinMax(1,1));

% Determine the expected value of the response for each ts
ta = ta_expectation3(ts,wm,N,dt,'wp',wp,'Type',estimator.type,'method','analytical','integrationMethod',method,'options',options);
T.ta = ta;
T.ts = ts;

% Find the bias
bias = mean((ta - ts).^2);

MinMax = [MinMax; min(ta)-3*wp*min(ta) max(ta)+3*wp*max(ta)];

% Find the variance
switch estimator.type
    case 'BLS'
        functionHandle = @(X,ts)integrand_BLS_quad(X,T,wm,wp,dt,N);
        vs = ndintegrate(functionHandle,MinMax,'method',method,'options',method_options,'ExtraVariables',T);
        variance = mean(vs);
        
    case 'sequential'
        % TODO
    case 'observer-actor'
        % TODO
        
end

%% Functions
function out = integrand_BLS_quad(X,ts,wm,wp,dt,N)
% Assign variables
tm = X(:,1:end-1);
tp = X(:,end);
method_opts.type = 'quad';
method_opts.dx = dt;
ta(:,1) = ta_expectation3(ts,wm,N,dt,'wp',wp,'Type',estimator.type,'method','analytical','integrationMethod','quad','options',options);
te = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts);

% Develop matrices
M = repmat(permute(tm,[2 3 1]),[1 1 1 size(tp,1) size(ts,1)]);
TE = repmat(permute(te,[2 3 1]),[1 1 1 size(tp,1) size(ts,1)]);
TP = repmat(permut(tp,[2 3 1]),[
TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
TP = repmat(tp,[1 size(ta,1)]);
TA = permute(ta',[size(tp,1) 1]);

% Find the probability functions
p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
p_tp_take_ta = ( (1./sqrt(2*pi)/wp/TE) .* exp( -(TP-TE).^2./(2*wp.^2.*TE.^2) ) );

% Determine the integrand for each combination of tm and tp
out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* p_tp_take_ta .* (TP - TA).^2;