function [tp, probTp, varargout] = ProbTp(ts,wm,wp,N,dt,varargin)
%% ta_expectation
%
%   ta = ta_expectation2(ts,wm,N)
%
%   Determines the expected value of the aim time across all possible
%   measurements based on the weber fraction of the measurement, wm, for N
%   measurements.
%
%   swe 20140616
%       Notes -
%           
%
%%

% Defaults
method_opts_default.dx = 0.01;
tpvec_default = 200:50:1400;
estimator_default.type = 'BLS';

% Parse inputs
p = inputParser;
addRequired(p,'ts');
addRequired(p,'wm');
addRequired(p,'wp');
addRequired(p,'N');
addRequired(p,'dt');
addParameter(p,'Type','BLS')
addParameter(p,'Support',NaN)
addParameter(p,'method','analytical');
addParameter(p,'method_options',method_opts_default);
addParameter(p,'trials',1000);
addParameter(p,'integrationMethod','quad');
addParameter(p,'options',NaN);
addParameter(p,'tpvec',tpvec_default);
addParameter(p,'estimator',estimator_default)
addParameter(p,'b',0)

parse(p,ts,wm,wp,N,dt,varargin{:})

ts = p.Results.ts;
wm = p.Results.wm;
wp = p.Results.wp;
N = p.Results.N;
dt = p.Results.dt;
Type = p.Results.Type;
Support = p.Results.Support;
method = p.Results.method;
method_options = p.Results.method_options;
trials = p.Results.trials;
integrationMethod = p.Results.integrationMethod;
options = p.Results.options;
tpvec = p.Results.tpvec;
estimator = p.Results.estimator;
b = p.Results.b;

if isnan(Support)
    tsmin = min(ts);
    tsmax = max(ts);
else
    tsmin = Support(1);
    tsmax = Support(2);
end

% Find the expected value of the estimate
switch method
    case 'analytical'
        
        estimator.wy = wp;
        probTp = p_tp(ts,tpvec(:),wm,wp,b,2,tsmin,tsmax,estimator);
        tp = NaN;

        
    case 'numerical'
        % Iterate for each ts
        %h = waitbar(0,['Simulating ' num2str(N) ' measurements']);
        tp = [];
        for i = 1:length(ts)
            % Generate measuments of each ts
            noise = wm*(ts(i)*ones(trials,N)).*randn(trials,N);
            tm = ts(i)*ones(trials,N) + noise;
            method_opts.type = 'quad';
            method_opts.dx = dt;
            BLS = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator);
            
            tp = [tp; BLS + wp*BLS.*randn(size(BLS))];
            
        end
        
        % Find the probability distribution
        N = hist(tp,tpvec);
        probTp = N/sum(N);
end
        

%% Functions
function p = p_tp(ts,tp,wm,wp,b,N,tsmin,tsmax,estimator)
%% For generating the probability of each tp
% [y1, y2] = ndgrid(ts,tp);
% Y = [y1(:) y2(:)];

Y{1} = ts;
Y{2} = tp;

functionHandle = @(tm,Y)(integrand(tm,Y,wm,wp,b,tsmin,tsmax,estimator));
options.dx = 0.001;
p = ndintegrate(functionHandle,repmat([tsmin-5*tsmin*wm tsmax+5*tsmax*wm],N,1),'method','quad','options',options,'ExtraVariables',Y);


function out = integrand(tm,Y,wm,wp,b,tsmin,tsmax,estimator)
%% For integration of across measurements in generating the p(tp|ts)

% Unpack extra variables
% ts = Y(:,1);
% tp = Y(:,2);
ts = Y{1};
tp = Y{2};

% Reshape the elements so that the computations can be performed correctly
% TS = repmat(permute(ts,[3 2 1]),size(tm,1),size(tm,2),1);
% TP = repmat(permute(tp,[3 2 1]),size(tm,1),1,1);
TS = repmat(ts,size(tm,1),1);
TP = repmat(permute(tp,[3 2 1]),size(tm,1),1,1);
method_opts.type = 'quad';
method_opts.dx = 10;
fBLS = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator) + b;
fBLS = repmat(fBLS,1,1,size(TP,3));

% Find the probability for each combination of variables
if size(tm,2) == 1
    p_tp_take_fBLS = (1./sqrt(2*pi*wp^2*fBLS.^2)).*exp(-(TP-fBLS).^2./(2*wp^2*fBLS.^2));
    p_tm_take_ts = (1./sqrt(2*pi*wm^2*TS.^2)).*exp(-sum((tm-TS).^2,2)./(2*wm^2*TS.^2));
elseif size(tm,2) == 2
%     TPmod = TP(:,:,1:end/2);
%     fBLSmod = fBLS(:,:,1:end/2);
%     p_tp_take_fBLS = (1./sqrt(2*pi*wp^2*fBLSmod.^2)).*exp(-(TPmod-fBLSmod).^2./(2*wp^2*fBLSmod.^2));
%     p_tp_take_fBLS = (1./sqrt(2*pi*wp^2*fBLS.^2)).*exp(-(TP-fBLS).^2./(2*wp^2*fBLS.^2));
%     TS1 = TS(:,1,1:2:end);
%     TS2 = TS(:,2,2:2:end);
%     p_tm_take_ts = (1./sqrt(2*pi*wm^2*TS1.^2)).*exp(-(repmat(tm(:,1),1,1,size(TS1,3))-TS1).^2./(2*wm^2*TS1.^2)) .* (1./sqrt(2*pi*wm^2*TS2.^2)).*exp(-(repmat(tm(:,2),1,1,size(TS2,3))-TS2).^2./(2*wm^2*TS2.^2)); %(1./(2*pi*wm^2*TS1.*TS2)).*exp(-(repmat(tm(:,1),1,1,size(TS1,3))-TS1).^2./(2*wm^2*TS1.^2) - (repmat(tm(:,2),1,1,size(TS2,3))-TS2).^2./(2*wm^2*TS2.^2));
%     p_tm_take_ts = (1./sqrt(2*pi*wm^2*TS(:,1).^2)).*exp(-(tm(:,1)-TS(:,1)).^2./(2*wm^2*TS(:,1).^2)) .* (1./sqrt(2*pi*wm^2*TS(:,2).^2)).*exp(-(tm(:,2)-TS(:,2)).^2./(2*wm^2*TS(:,2).^2));
    wm2 = wm;
    if isfield(estimator,'wm_drift')
        wm1 = estimator.wm_drift;
    else
        wm1 = wm;
    end
%     TPmod = TP(:,:,1:end/2);
%     fBLSmod = fBLS(:,:,1:end/2);
%     p_tp_take_fBLS = (1./sqrt(2*pi*wp^2*fBLSmod.^2)).*exp(-(TPmod-fBLSmod).^2./(2*wp^2*fBLSmod.^2));
    p_tp_take_fBLS = (1./sqrt(2*pi*wp^2*fBLS.^2)).*exp(-(TP-fBLS).^2./(2*wp^2*fBLS.^2));
%     TS1 = TS(:,1,1:2:end);
%     TS2 = TS(:,2,2:2:end);
%     p_tm_take_ts = (1./sqrt(2*pi*wm^2*TS1.^2)).*exp(-(repmat(tm(:,1),1,1,size(TS1,3))-TS1).^2./(2*wm^2*TS1.^2)) .* (1./sqrt(2*pi*wm^2*TS2.^2)).*exp(-(repmat(tm(:,2),1,1,size(TS2,3))-TS2).^2./(2*wm^2*TS2.^2)); %(1./(2*pi*wm^2*TS1.*TS2)).*exp(-(repmat(tm(:,1),1,1,size(TS1,3))-TS1).^2./(2*wm^2*TS1.^2) - (repmat(tm(:,2),1,1,size(TS2,3))-TS2).^2./(2*wm^2*TS2.^2));
    p_tm_take_ts = (1./sqrt(2*pi*wm1^2*TS(:,1).^2)).*exp(-(tm(:,1)-TS(:,1)).^2./(2*wm1^2*TS(:,1).^2)) .* ...
        (1./sqrt(2*pi*wm2^2*TS(:,2).^2)).*exp(-(tm(:,2)-TS(:,2)).^2./(2*wm2^2*TS(:,2).^2));
end

out = permute(p_tp_take_fBLS.*repmat(p_tm_take_ts,1,1,size(TP,3)),[1 3 2]);

%% Functions
% BLS
% Simpson's quad
% function out = integrandFunction_quad(tm,ts,tsmin,tsmax,wm,dt,N)
% 
% method_opts.type = 'quad';
% method_opts.dx = dt;
% TA = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts);
% 
% M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
% TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
% %p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
% p_tm_take_ts = ( (1./sqrt(2*pi)./wm./TS(1,:,:,:)).^N .* exp( -(sum((TS-M).^2,1))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
% %p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
% out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA,[1 size(ts)]);
% 
% % Monte Carlo integration
% function out = integrandFunction_MonteCarlo(tm,ts,tsmin,tsmax,wm,Npoints,N)
% 
% method_opts.type = 'MonteCarlo';
% method_opts.N = Npoints;
% TA = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts);
% 
% M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
% TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
% %p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
% p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(sum((TS-M).^2,1))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
% %p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
% out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA',[1 size(ts)]);
% 
% % Monte Carlo integration, in batch
% function out = integrandFunction_MonteCarlo_batch(tm,ts,tsmin,tsmax,wm,Npoints,batch_sz,N)
% 
% method_opts.type = 'MonteCarlo_batch';
% method_opts.N = Npoints;
% method_opts.batch_sz = batch_sz;
% TA = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts);
% 
% M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
% TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
% %p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
% p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(sum((TS-M).^2,1))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
% %p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
% out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA',[1 size(ts)]);
% 
% % Sequential
% % Simpson's quad
% function out = integrandFunction_quadSequential(tm,ts,tsmin,tsmax,wm,dt,N,deterministicModel,probabilisticModel,boundry)
% method_opts.type = 'quad';
% method_opts.dx = dt;
% estimator_opts.type = 'sequential';
% estimator_opts.deterministicModel = deterministicModel;
% estimator_opts.probabilisticModel = probabilisticModel;
% estimator_opts.boundry = boundry;
% %estimator_opts.p = params;
% TA = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator_opts);
% 
% M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
% TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
% %p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
% p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(sum((TS-M).^2,1))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
% %p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
% out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA,[1 size(ts)]);
% 
% % gBLS
% % Simpson's quad
% function out = integrandFunction_gBLS_quad(tm,ts,tsmin,tsmax,wm,g,dt,N)
% 
% method_opts.type = 'quad';
% method_opts.dx = dt;
% TA = g*ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts);
% 
% M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
% TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
% %p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
% p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(sum((TS-M).^2,1))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
% %p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
% out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA,[1 size(ts)]);