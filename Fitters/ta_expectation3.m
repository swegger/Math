function [ta, ta_std, varargout] = ta_expectation3(ts,wm,N,dt,varargin)
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

% Parse inputs
p = inputParser;
addRequired(p,'ts');
addRequired(p,'wm');
addRequired(p,'N');
addRequired(p,'dt');
addParameter(p,'Type','BLS')
addParameter(p,'wp',NaN);
addParameter(p,'Support',NaN)
addParameter(p,'method','analytical');
addParameter(p,'method_options',method_opts_default);
addParameter(p,'trials',1000);
addParameter(p,'integrationMethod','quad');
addParameter(p,'options',NaN);
addParameter(p,'sigp',0);


parse(p,ts,wm,N,dt,varargin{:})

ts = p.Results.ts;
wm = p.Results.wm;
N = p.Results.N;
dt = p.Results.dt;
Type = p.Results.Type;
wp = p.Results.wp;
Support = p.Results.Support;
method = p.Results.method;
method_options = p.Results.method_options;
trials = p.Results.trials;
integrationMethod = p.Results.integrationMethod;
options = p.Results.options;
sigp = p.Results.sigp;

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
        switch Type
            case 'BLS'
                switch integrationMethod
                    case 'integral'
                        error('integrationMethod "integral" not yet supported!')
                    case 'quad'
                        functionHandle = @(tm,ts)(integrandFunction_quad(tm,ts,tsmin,tsmax,wm,dt,N));
                        %method_options.dx = 0.01;
                        ta = ndintegrate(functionHandle,repmat([tsmin-2*wm*tsmin tsmax+3*wm*tsmax],N,1),'method','quad','options',method_options,'ExtraVariables',ts);
                        
                    case 'MonteCarlo'
                        functionHandle = @(tm,ts)(integrandFunction_MonteCarlo(tm,ts,tsmin,tsmax,wm,method_options.N,N));
                        %method_options.N = 100000;
                        ta = ndintegrate(functionHandle,repmat([tsmin-3*wm*tsmin tsmax+3*wm*tsmax],N,1),'method','MonteCarlo','options',method_options,'ExtraVariables',ts);
                        %error('integrationMethod "MonteCarlo" not yet supported!')
                        
                    case 'MonteCarlo_batch'
                        functionHandle = @(tm,ts)(integrandFunction_MonteCarlo_batch(tm,ts,tsmin,tsmax,wm,method_options.N,method_options.batch_sz,N));
                        %method_options.N = 100000;
                        ta = ndintegrate(functionHandle,repmat([tsmin-3*wm*tsmin tsmax+3*wm*tsmax],N,1),'method','MonteCarlo_batch','options',method_options,'ExtraVariables',ts);
                        
                    otherwise
                        error(['integrationMethod ' integrationMethod ' not recognized!'])
                        
                end
                
            case 'MAP'
                switch integrationMethod
                    case 'integral'
                        error('integrationMethod "integral" not yet supported!')
                    case 'quad'
                        functionHandle = @(tm,ts)(MAP_integrandFunction_quad(tm,ts,tsmin,tsmax,wm,dt,N));
                        %method_options.dx = 0.01;
                        ta = ndintegrate(functionHandle,repmat([tsmin-2*wm*tsmin tsmax+3*wm*tsmax],N,1),'method','quad','options',method_options,'ExtraVariables',ts);
                        
                    case 'MonteCarlo'
                        error(['integrationMethod ' integrationMethod ' not yet supported!'])
                        functionHandle = @(tm,ts)(integrandFunction_MonteCarlo(tm,ts,tsmin,tsmax,wm,method_options.N,N));
                        %method_options.N = 100000;
                        ta = ndintegrate(functionHandle,repmat([tsmin-3*wm*tsmin tsmax+3*wm*tsmax],N,1),'method','MonteCarlo','options',method_options,'ExtraVariables',ts);
                        %error('integrationMethod "MonteCarlo" not yet supported!')
                        
                    case 'MonteCarlo_batch'
                        error(['integrationMethod ' integrationMethod ' not yet supported!'])
                        functionHandle = @(tm,ts)(integrandFunction_MonteCarlo_batch(tm,ts,tsmin,tsmax,wm,method_options.N,method_options.batch_sz,N));
                        %method_options.N = 100000;
                        ta = ndintegrate(functionHandle,repmat([tsmin-3*wm*tsmin tsmax+3*wm*tsmax],N,1),'method','MonteCarlo_batch','options',method_options,'ExtraVariables',ts);
                        
                    otherwise
                        error(['integrationMethod ' integrationMethod ' not recognized!'])
                        
                end
                
                
            case 'ObsAct'
                switch integrationMethod
                    case 'integral'
                        error('integrationMethod "integral" not yet supported!')
                    case 'quad'
                        functionHandle = @(tm,ts)(ObsAct_integrandFunction_quad(tm,ts,tsmin,tsmax,wm,dt,N,wp));
                        %method_options.dx = 0.01;
                        ta = ndintegrate(functionHandle,repmat([tsmin-2*wm*tsmin tsmax+3*wm*tsmax],N,1),'method','quad','options',method_options,'ExtraVariables',ts);
                        
                    case 'MonteCarlo'
                        functionHandle = @(tm,ts)(ObsAct_integrandFunction_MonteCarlo(tm,ts,tsmin,tsmax,wm,method_options.N,N,wp));
                        %method_options.N = 100000;
                        ta = ndintegrate(functionHandle,repmat([tsmin-3*wm*tsmin tsmax+3*wm*tsmax],N,1),'method','MonteCarlo','options',method_options,'ExtraVariables',ts);
                        %error('integrationMethod "MonteCarlo" not yet supported!')
                        
                    case 'MonteCarlo_batch'
                        functionHandle = @(tm,ts)(ObsAct_integrandFunction_MonteCarlo_batch(tm,ts,tsmin,tsmax,wm,method_options.N,method_options.batch_sz,N,wp));
                        %method_options.N = 100000;
                        ta = ndintegrate(functionHandle,repmat([tsmin-3*wm*tsmin tsmax+3*wm*tsmax],N,1),'method','MonteCarlo_batch','options',method_options,'ExtraVariables',ts);
                        
                    otherwise
                        error(['integrationMethod ' integrationMethod ' not recognized!'])
                        
                end
                
                
            case 'gBLS'
                switch integrationMethod
                    case 'integral'
                        error('integrationMethod "integral" not yet supported!')
                    case 'quad'
                        functionHandle = @(tm,ts)(integrandFunction_gBLS_quad(tm,ts,tsmin,tsmax,wm,options.g,dt,N));
                        %method_options.dx = 0.01;
                        ta = ndintegrate(functionHandle,repmat([tsmin-3*wm*tsmin tsmax+3*wm*tsmax],N,1),'method','quad','options',method_options,'ExtraVariables',ts);
                        
                    case 'MonteCarlo'
                        error('integrationMethod "MonteCarlo" not yet supported!')
                end
                
            case 'sequential'
                switch integrationMethod
                    case 'integral'
                        error('integrationMethod "integral" not yet supported!')
                    case 'quad'
                        deterministicModel = options.deterministicModel;
                        probabilisticModel = options.probabilisticModel;
                        boundry = options.boundry;
                        %params = options.parameters;
                        functionHandle = @(tm,ts)(integrandFunction_quadSequential(tm,ts,tsmin,tsmax,wm,dt,N,deterministicModel,probabilisticModel,boundry));
                        %method_options.dx = 0.01;
                        ta = ndintegrate(functionHandle,repmat([tsmin-3*wm*tsmin tsmax+3*wm*tsmax],N,1),'method','quad','options',method_options,'ExtraVariables',ts);
                    case 'trapz'
                        error('integrationMethod "trapz" not yet supported!')
                    case 'MonteCarlo'
                        error('integrationMethod "MonteCarlo" not yet supported!')
                end
                
            case 'aveMeasurements'
                error('Analytical solution to aveMeasurement model not yet supported!')
                
            otherwise
                error(['Analytical solution to ' Type ' model not yet supported'])
        end
        
    case 'numerical'
        switch Type
            case 'BLS'
                % Iterate for each ts
                %h = waitbar(0,['Simulating ' num2str(N) ' measurements']);
                errors = zeros(trials,length(ts));
                for i = 1:length(ts)
                    % Generate measuments of each ts
                    noise = wm*(ts(i)*ones(trials,N)).*randn(trials,N);
                    tm = ts(i)*ones(trials,N) + noise;
                    method_opts.type = 'quad';
                    method_opts.dx = dt;
                    BLS = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts);
                    
                    BLS = BLS + wp*BLS.*randn(size(BLS));
                    
                    errors(:,i) = BLS - ts(i);
                    ta(i) = mean(BLS);
                    ta_std(i) = std(BLS);
                    
                    %waitbar(i/length(ts))
                end
                
                % Bias
                rmse = sqrt(mean(errors(:).^2));
                bias2 = mean((ta(:)-ts(:)).^2);
                v = mean(ta_std.^2);
                varargout{1} = bias2;
                varargout{2} = v;
                varargout{3} = rmse;
            
            case 'MAP'
                % Iterate for each ts
                %h = waitbar(0,['Simulating ' num2str(N) ' measurements']);
                errors = zeros(trials,length(ts));
                for i = 1:length(ts)
                    % Generate measuments of each ts
                    noise = wm*(ts(i)*ones(trials,N)).*randn(trials,N);
                    tm = ts(i)*ones(trials,N) + noise;
                    method_opts.type = 'analytical';
                    method_opts.dx = dt;
                    estimator.type = 'MAP';
                    MAP = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator);
                    
                    MAP = MAP + wp*MAP.*randn(size(MAP));
                    
                    errors(:,i) = MAP - ts(i);
                    ta(i) = mean(MAP);
                    ta_std(i) = std(MAP);
                    
                    %waitbar(i/length(ts))
                end
                
                % Bias
                rmse = sqrt(mean(errors(:).^2));
                bias2 = mean((ta(:)-ts(:)).^2);
                v = mean(ta_std.^2);
                varargout{1} = bias2;
                varargout{2} = v;
                varargout{3} = rmse;
                
            case 'ObsAct'
                % Iterate for each ts
                %h = waitbar(0,['Simulating ' num2str(N) ' measurements']);
                errors = zeros(trials,length(ts));
                for i = 1:length(ts)
                    % Generate measuments of each ts
                    noise = wm*(ts(i)*ones(trials,N)).*randn(trials,N);
                    tm = ts(i)*ones(trials,N) + noise;
                    method_opts.type = 'quad';
                    method_opts.dx = dt;
                    estimator.type = 'ObsAct';
                    estimator.wy = wp;
                    BLS = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator);
                    
                    BLS = BLS + wp*BLS.*randn(size(BLS));
                    
                    errors(:,i) = BLS - ts(i);
                    ta(i) = mean(BLS);
                    ta_std(i) = std(BLS);
                    
                    %waitbar(i/length(ts))
                end
                
                % Bias
                rmse = sqrt(mean(errors(:).^2));
                bias2 = mean((ta(:)-ts(:)).^2);
                v = mean(ta_std.^2);
                varargout{1} = bias2;
                varargout{2} = v;
                varargout{3} = rmse;                
                
            case 'gBLS'
                error('gBLS numerical simulation not yet supported!')
                
            case 'sequential'
                error('sequential numerical simulation not yet supported!')
                
            case 'aveMeasurements'
                
                % Iterate for each ts
                %h = waitbar(0,['Simulating ' num2str(N) ' measurements']);
                errors = zeros(trials,length(ts));
                for i = 1:length(ts)
                    % Generate measuments of each ts
                    noise = wm*(ts(i)*ones(trials,N)).*randn(trials,N);
                    tm = ts(i)*ones(trials,N) + noise;
                    method_opts.type = 'quad';
                    method_opts.dx = dt;
                    estimator.type = 'weightedMean';
                    estimator.weights = ones(1:N)/N;
                    E = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator);
                    
                    E = E + wp*E.*randn(size(E));
                    
                    errors(:,i) = E - ts(i);
                    ta(i) = mean(E);
                    ta_std(i) = std(E);
                    
                    %waitbar(i/length(ts))
                end
                
                % Bias
                rmse = sqrt(mean(errors(:).^2));
                bias2 = mean((ta(:)-ts(:)).^2);
                v = mean(ta_std.^2);
                varargout{1} = bias2;
                varargout{2} = v;
                varargout{3} = rmse;
                
            case 'aveEstimates'
                
                % Iterate for each ts
                %h = waitbar(0,['Simulating ' num2str(N) ' measurements']);
                errors = zeros(trials,length(ts));
                for i = 1:length(ts)
                    % Generate measuments of each ts
                    noise = wm*(ts(i)*ones(trials,N)).*randn(trials,N);
                    tm = ts(i)*ones(trials,N) + noise;
                    method_opts.type = 'quad';
                    method_opts.dx = dt;
                    estimator.type = 'aveEstimates';
                    estimator.weights = ones(1:N)/N;
                    E = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator);
                    
                    E = E + wp*E.*randn(size(E));
                    
                    errors(:,i) = E - ts(i);
                    ta(i) = mean(E);
                    ta_std(i) = std(E);
                    
                    %waitbar(i/length(ts))
                end
                
                % Bias
                rmse = sqrt(mean(errors(:).^2));
                bias2 = mean((ta(:)-ts(:)).^2);
                v = mean(ta_std.^2);
                varargout{1} = bias2;
                varargout{2} = v;
                varargout{3} = rmse;
                
            case 'BLS_wm_wp_sigp'
                % Iterate for each ts
                %h = waitbar(0,['Simulating ' num2str(N) ' measurements']);
                errors = zeros(trials,length(ts));
                for i = 1:length(ts)
                    % Generate measuments of each ts
                    noise = wm*(ts(i)*ones(trials,N)).*randn(trials,N);
                    tm = ts(i)*ones(trials,N) + noise;
                    method_opts.type = 'quad';
                    method_opts.dx = dt;
                    BLS = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts);
                    
                    BLS = BLS + wp*BLS.*randn(size(BLS)) + sigp*randn(size(BLS));
                    
                    errors(:,i) = BLS - ts(i);
                    ta(i) = mean(BLS);
                    ta_std(i) = std(BLS);
                    
                    %waitbar(i/length(ts))
                end
                
                % Bias
                rmse = sqrt(mean(errors(:).^2));
                bias2 = mean((ta(:)-ts(:)).^2);
                v = mean(ta_std.^2);
                varargout{1} = bias2;
                varargout{2} = v;
                varargout{3} = rmse;
                
            otherwise
                error(['Numerical simulation of model type ' Type ' not yet supported!'])
        end
end
        

%% Functions
% BLS
% Simpson's quad
function out = integrandFunction_quad(tm,ts,tsmin,tsmax,wm,dt,N)

method_opts.type = 'quad';
method_opts.dx = dt;
TA = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts);

M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
%p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
p_tm_take_ts = ( (1./sqrt(2*pi)./wm./TS(1,:,:,:)).^N .* exp( -(sum((TS-M).^2,1))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
%p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA,[1 size(ts)]);

% Monte Carlo integration
function out = integrandFunction_MonteCarlo(tm,ts,tsmin,tsmax,wm,Npoints,N)

method_opts.type = 'MonteCarlo';
method_opts.N = Npoints;
TA = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts);

M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
%p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(sum((TS-M).^2,1))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
%p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA',[1 size(ts)]);

% Monte Carlo integration, in batch
function out = integrandFunction_MonteCarlo_batch(tm,ts,tsmin,tsmax,wm,Npoints,batch_sz,N)

method_opts.type = 'MonteCarlo_batch';
method_opts.N = Npoints;
method_opts.batch_sz = batch_sz;
TA = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts);

M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
%p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(sum((TS-M).^2,1))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
%p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA',[1 size(ts)]);

% MAP
% Simpson's quad
function out = MAP_integrandFunction_quad(tm,ts,tsmin,tsmax,wm,dt,N)

method_opts.type = 'analytical';
method_opts.dx = dt;
estimator.type = 'MAP';
TA = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator);

M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
%p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
p_tm_take_ts = ( (1./sqrt(2*pi)./wm./TS(1,:,:,:)).^N .* exp( -(sum((TS-M).^2,1))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
%p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA,[1 size(ts)]);

% ObsAct
% Simpson's quad
function out = ObsAct_integrandFunction_quad(tm,ts,tsmin,tsmax,wm,dt,N,wy)

method_opts.type = 'quad';
method_opts.dx = dt;
estimator.type = 'ObsAct';
estimator.wy = wy;
TA = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator);

M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
%p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
p_tm_take_ts = ( (1./sqrt(2*pi)./wm./TS(1,:,:,:)).^N .* exp( -(sum((TS-M).^2,1))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
%p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA,[1 size(ts)]);

% Monte Carlo integration
function out = ObsAct_integrandFunction_MonteCarlo(tm,ts,tsmin,tsmax,wm,Npoints,N,wy)

method_opts.type = 'MonteCarlo';
method_opts.N = Npoints;
estimator.type = 'ObsAct';
estimator.wy = wy;
TA = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator);

M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
%p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(sum((TS-M).^2,1))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
%p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA',[1 size(ts)]);

% Monte Carlo integration, in batch
function out = ObsAct_integrandFunction_MonteCarlo_batch(tm,ts,tsmin,tsmax,wm,Npoints,batch_sz,N,wy)

method_opts.type = 'MonteCarlo_batch';
method_opts.N = Npoints;
method_opts.batch_sz = batch_sz;
estimator.type = 'ObsAct';
estimator.wy = wy;
TA = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator);

M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
%p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(sum((TS-M).^2,1))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
%p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA',[1 size(ts)]);

% Sequential
% Simpson's quad
function out = integrandFunction_quadSequential(tm,ts,tsmin,tsmax,wm,dt,N,deterministicModel,probabilisticModel,boundry)
method_opts.type = 'quad';
method_opts.dx = dt;
estimator_opts.type = 'sequential';
estimator_opts.deterministicModel = deterministicModel;
estimator_opts.probabilisticModel = probabilisticModel;
estimator_opts.boundry = boundry;
%estimator_opts.p = params;
TA = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator_opts);

M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
%p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(sum((TS-M).^2,1))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
%p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA,[1 size(ts)]);

% gBLS
% Simpson's quad
function out = integrandFunction_gBLS_quad(tm,ts,tsmin,tsmax,wm,g,dt,N)

method_opts.type = 'quad';
method_opts.dx = dt;
TA = g*ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts);

M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
%p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(sum((TS-M).^2,1))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
%p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA,[1 size(ts)]);