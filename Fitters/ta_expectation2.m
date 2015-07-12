function [ta ta_std varargout] = ta_expectation2(ts,wm,N,varargin)
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

% Parse inputs
p = inputParser;
addRequired(p,'ts');
addRequired(p,'wm');
addRequired(p,'N');
addParamValue(p,'Type','BLS')
addParamValue(p,'wp',NaN);
addParamValue(p,'method','analytical');
addParamValue(p,'trials',1000);
addParamValue(p,'integrationMethod','quad');
addParamValue(p,'options',NaN);


parse(p,ts,wm,N,varargin{:})

ts = p.Results.ts;
wm = p.Results.wm;
N = p.Results.N;
Type = p.Results.Type;
wp = p.Results.wp;
method = p.Results.method;
trials = p.Results.trials;
integrationMethod = p.Results.integrationMethod;
options = p.Results.options;

% 
% 
% switch method
%     case 'analytical'
%         if N == 1
%             switch Type
%                 case 'BLS'
%                     % BLS solution for a given tm
%                     f = @(m,xmin,xmax,wm)(integral(@(x)(x.*exp(-1./2.*(x-m).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)(exp(-1./2.*(x-m).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true));
%                     
%                     % Find the average ta for each ts by integrating f(tm)*p(tm|ts) over all tm
%                     integrand = @(ts,tm,wm,tsmin,tsmax)( (1./sqrt(2*pi*wm.^2*ts.^2)) .* exp(-1./2.*(tm-ts).^2./(wm.*ts).^2) .* f(tm,tsmin,tsmax,wm) );
%                     ta = integral(@(tm)(integrand(ts,tm,wm(1),ts(1),ts(end))),ts(1)-15*wm(1)*ts(1),ts(end)+15*wm(1)*ts(end),'ArrayValued',true);
%                     
%                 case 'ObsAct'
%                     if isnan(wp)
%                         error('wp not supplied')
%                     else
%                         % BLS solution for a given tm
%                         f = @(m,xmin,xmax,wm)( integral(@(x)(x.*exp(-1./2.*(x-m).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)(exp(-1./2.*(x-m).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true) ./ (1+wp.^2) );
%                         
%                         % Find the average ta for each ts by integrating f(tm)*p(tm|ts) over all tm
%                         integrand = @(ts,tm,wm,tsmin,tsmax)( (1./sqrt(2*pi*wm.^2*ts.^2)) .* exp(-1./2.*(tm-ts).^2./(wm.*ts).^2) .* f(tm,tsmin,tsmax,wm) );
%                         ta = integral(@(tm)(integrand(ts,tm,wm(1),ts(1),ts(end))),ts(1)-15*wm(1)*ts(1),ts(end)+15*wm(1)*ts(end),'ArrayValued',true);
%                     end
%             end
%                     
%             
%         elseif N == 2
%             switch Type
%                 case 'BLS'
%                     % BLS solution for a given combination of tm1 and tm2
%                     f = @(m1,m2,xmin,xmax,wm)(integral(@(x)(x.*(1./sqrt(2.*pi.*wm.^2.*x.^2)).^2.*exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)((1./sqrt(2.*pi.*wm.^2.*x.^2)).^2.*exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true));
%                     
%                     % Find the average ta for each ts by integrating f(tm1,tm2)*p(tm1,tm2|ts) over all tm1 and tm2
%                     integrand = @(ts,tm1,tm2,wm,tsmin,tsmax)( (1./sqrt(2.*pi.*wm.^2.*ts.^2)).^2.*exp(-1./2.*(ts-tm1).^2./wm.^2./ts.^2).*exp(-1./2.*(ts-tm2).^2./wm.^2./ts.^2) .* f(tm1,tm2,tsmin,tsmax,wm) );
%                     h = waitbar(0,[num2str(N) ' measurements']);
%                     for i = 1:length(ts)
%                         tsi = ts(i);
%                         ta(i) = integral2(@(tm1,tm2)(integrand(tsi,tm1,tm2,wm(1),ts(1),ts(end))),ts(1)-15*wm(1)*ts(1),ts(end)+15*wm(1)*ts(end),ts(1)-15*wm(1)*ts(1),ts(end)+15*wm(1)*ts(end));
%                         waitbar(i/length(ts))
%                     end
%                     close(h)
%                 
%                 case 'ObsAct'
%                     if isnan(wp)
%                         error('wp not supplied')
%                     else
%                         % BLS solution for a given combination of tm1 and tm2
%                         f = @(m1,m2,xmin,xmax,wm)( integral(@(x)(x.*(1./sqrt(2.*pi.*wm.^2.*x.^2)).^2.*exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)((1./sqrt(2.*pi.*wm.^2.*x.^2)).^2.*exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true) ./ (1+wp.^2) );
%                         
%                         % Find the average ta for each ts by integrating f(tm1,tm2)*p(tm1,tm2|ts) over all tm1 and tm2
%                         integrand = @(ts,tm1,tm2,wm,tsmin,tsmax)( (1./sqrt(2.*pi.*wm.^2.*ts.^2)).^2.*exp(-1./2.*(ts-tm1).^2./wm.^2./ts.^2).*exp(-1./2.*(ts-tm2).^2./wm.^2./ts.^2) .* f(tm1,tm2,tsmin,tsmax,wm) );
%                         h = waitbar(0,[num2str(N) ' measurements']);
%                         for i = 1:length(ts)
%                             tsi = ts(i);
%                             ta(i) = integral2(@(tm1,tm2)(integrand(tsi,tm1,tm2,wm(1),ts(1),ts(end))),ts(1)-15*wm(1)*ts(1),ts(end)+15*wm(1)*ts(end),ts(1)-15*wm(1)*ts(1),ts(end)+15*wm(1)*ts(end));
%                             waitbar(i/length(ts))
%                         end
%                         close(h)
%                     end
%                     
%             end
%             
%         elseif N == 3
%             % Analytical solultion: TODO
%             % BLS solution for a given combination of tm1-3
%             f = @(V,xmin,xmax,wm)(integral(@(x)(x.*(1./(sqrt(2*pi)*wm*x)).^3 .* exp( -V ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true)./integral(@(x)((1./(sqrt(2*pi)*wm*x)).^3 .* exp( -V ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true));
%             
%             % Find the average ta for each ts by numerical intgraiton of f(tm_vec)*p(tm_vec|ts) over all tm_vec
%             integrand = @(ts,V,wm,tsmin,tsmax)( ( (1./(sqrt(2*pi)*wm*ts)).^3 .* exp( -V ./ (2*wm.^2.*ts.^2) ) ) .* f(V,tsmin,tsmax,wm) );
%             tm_vec = ts(1)-100:10:ts(end)+100;
%             [TM(1,1,:,:,:) TM(2,1,:,:,:) TM(3,1,:,:,:)] = ndgrid(tm_vec);
%             
%             for i = 1:length(ts)
%                 tic
%                 V = mmx('mult',TM-ts(i),TM-ts(i),'t');
%                 TA = integrand(ts(i),V,wm,ts(1),ts(end));
%                 ta(i) = trapz(trapz(trapz(TA)));
%                 toc
%             end
%             
% %             for i = 1:length(ts)
% %                 tsi = ts(i);
% %                 TA = NaN(length(tm_vec),length(tm_vec),length(tm_vec));
% %                 
% %                 for ii = 1:length(tm_vec)
% %                     for jj = 1:length(tm_vec)
% %                         tic
% %                         for kk = 1:length(tm_vec)
% %                             TA(ii,jj,kk) = integrand(tsi,V,wm,ts(1),ts(end));
% %                         end
% %                         toc
% %                     end
% %                 end
% %                 ta(i) = trapz(trapz(trapz(TA)));
% %             end
%             
%             %error('Analytical solution for N == 3 measurements not supported.')
%             
%         else
%             % TODO
%             error('Analytical solution for N >= 3 measurements not supported.')
%             
%         end
%         
%         varargout{1} = NaN;
%         varargout{2} = NaN;
% 
%     case 'numerical'
%         % Iterate for each ts
%         h = waitbar(0,['Simulating ' num2str(N) ' measurements']);
%         for i = 1:length(ts)
%             % Generate measuments of each ts
%             noise = wm*(ts(i)*ones(trials,N)).*randn(trials,N);
%             tm = ts(i)*ones(trials,N) + noise;
%             
%             % Create expression for the likelihood
%             Likelihoodf = @(w,N,tm,ts)( (1./(sqrt(2*pi)*w*ts)).^N .* exp( -(tm-ts)'*(tm-ts) ./ (2*w.^2.*ts.^2) ) );
%             
%             % Use Bayes' theorem to find the BLS solution for each trial
%             % Find the normalization factor
%             for j = 1:trials
%                 
%                 for k = 1:length(ts)
%                     P_tm_take_ts(k) = Likelihoodf(wm,N,tm(j,:)',ts(k));
%                 end
%                 Z = trapz(P_tm_take_ts);
%                 
%                 % Now calculate the BLS solution
%                 for k = 1:length(ts)
%                     L(k) = Likelihoodf(wm,N,tm(j,:)',ts(k));
%                 end
%                 BLS(j,:) = trapz(ts.*L/Z);
%                 if ~isnan(wp)
%                     BLS(j,:) = BLS(j,:) + wp*BLS(j,:)*randn;
%                 end
%             end
%             ta(i) = mean(BLS);
%             ta_std(i) = std(BLS);
%             
%             waitbar(i/length(ts))
%         end
%         
%         % Bias
%         bias2 = mean((ta-ts).^2);
%         v = mean(ta_std.^2);
%         varargout{1} = bias2;
%         varargout{2} = v;
%         
%         close(h)
% end


switch method
    case 'analytical'
        switch Type
            case 'BLS'
                switch integrationMethod
                    case 'integral'
                        error('integrationMethod "integral" not yet supported!')
                    case 'quad'
                        functionHandle = @(tm,ts)(integrandFunction_quad(tm,ts,min(ts),max(ts),wm,1,N));
                        method_options.dx = 0.01;
                        ta = ndintegrate(functionHandle,repmat([100 2500],N,1),'method','quad','options',method_options,'ExtraVariables',ts);
                        
                    case 'MonteCarlo'
                        functionHandle = @(tm,ts)(integrandFunction_MonteCarlo(tm,ts,min(ts),max(ts),wm,1,N));
                        method_options.N = 100000;
                        ta = ndintegrate(functionHandle,repmat([100 2500],N,1),'method','MonteCarlo','options',method_options,'ExtraVariables',ts);
                        %error('integrationMethod "MonteCarlo" not yet supported!')
                end
                
            case 'sequential'
                switch integrationMethod
                    case 'integral'
                        error('integrationMethod "integral" not yet supported!')
                    case 'quad'
                        transitionModel = options.transitionModel;
                        params = options.parameters;
                        functionHandle = @(tm,ts)(integrandFunction_quadSequential(tm,ts,min(ts),max(ts),wm,1,N,transitionModel,params));
                        method_options.dx = 0.01;
                        ta = ndintegrate(functionHandle,repmat([100 2500],N,1),'method','quad','options',method_options,'ExtraVariables',ts);
                    case 'trapz'
                        error('integrationMethod "trapz" not yet supported!')
                    case 'MonteCarlo'
                        error('integrationMethod "MonteCarlo" not yet supported!')
                end
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
p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA,[1 size(ts)]);

% Monte Carlo integration
function out = integrandFunction_MonteCarlo(tm,ts,tsmin,tsmax,wm,Npoints,N)

method_opts.type = 'MonteCarlo';
method_opts.N = Npoints;
TA = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts);

M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
%p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA,[1 size(ts)]);

% Sequential
% Simpson's quad
function out = integrandFunction_quadSequential(tm,ts,tsmin,tsmax,wm,dt,N,transitionModel,params)
method_opts.type = 'quad';
method_opts.dx = dt;
estimator_opts.type = 'sequential';
estimator_opts.transitionModel = transitionModel;
estimator_opts.p = params;
TA = ScalarBayesEstimators(tm,wm,tsmin,tsmax,'method',method_opts,'estimator',estimator_opts);

M = repmat(permute(tm,[2 3 1]),[1 1 1 size(ts,1)]);
TS = repmat(permute(ts,[4 3 2 1]),[size(M,1) size(M,2) size(M,3) 1]);
%p_tm_take_ts = ( (1./sqrt(2*pi)*wm*X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) )
p_tm_take_ts = ( (1./sqrt(2*pi)/wm/TS(1,:,:,:)).^N .* exp( -(mmx('mult',permute(TS-M,[2 1 3 4]),TS-M))./(2*wm.^2.*TS(1,:,:,:).^2) ) );
out = squeeze(permute(p_tm_take_ts,[3 4 1 2])) .* repmat(TA',[1 size(ts)]);
