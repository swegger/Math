function [wm, wy, llikelihood] = ObsAct_fitter(x,y,varargin)
%% FITBAYESOBSERVERMODEL
%
%   Fits the Baysian Observer Model (Jazayeri and Shadlen, 2010) with a bias 
%   to time of reproduction and sample time data.
%
%   Model is in the form:
%       p(y|x,wm,wy) = Integral(p(y|f(m),wy)p(m|x,wm)dm)
%       p(m|x,wm) = (1/sqrt(2*pi*(wm*x)^2))*exp(-(x-m)^2/(2*(wm*x)^2)
%       p(y|f(m),wy) = (1/sqrt(2*pi*(wy*f(m))^2))*exp(-(y-f(m))^2/(2*(wy*f(m)^2)
%
%
%%

% Parse input
ICflg = 0;
FitterFlg = 0;
for i = 1:length(varargin)
    if isstr(varargin{i})
        if strcmp(varargin{i},'InitCond')
            ICflg = 1;
            ICnum = i;
        elseif strcmp(varargin{i},'FitType')
            FitterFlg = 1;
            Fitnum = i;
        end
    end
end

if ICflg
    IC = varargin{ICnum+1};
else
    IC = [0.1 0.06 0];
end

if FitterFlg
    FitType = varargin{Fitnum+1};
    switch FitType
        case 'integral'
            xmin = varargin{Fitnum+2}(1);
            xmax = varargin{Fitnum+2}(2);
        
        case 'trapz'
            m = varargin{Fitnum+2};
    end
        
else
    FitType = 'integral';
    xmin = min(x);
    xmax = max(x);
end

if nargin < 3
    error('Not enought input arguments!')
end


% Fit according to minimizer
OPTIONS = optimset('Display','iter');

switch FitType
    case 'integral'
        % Use integral
        BLSminimizant = @(p)logLikelihood(p(1),p(2),x,y,xmin,xmax);
        
    case 'trapz'
        % Use trapz
        BLSminimizant = @(p)logLikelihoodTRAPZ(p(1),p(2),m,x,y);
end

minimizer = 'fminsearch(BLSminimizant, [wM_ini wP_ini], OPTIONS);';
wM_ini = IC(1);
wP_ini = IC(2);
[lparams llikelihood] = eval(minimizer);

wm = lparams(1);
wy = lparams(2);
llikelihood = llikelihood;
    

%% Function to be minimized
function logL = logLikelihood(wm,wy,x,y,xmin,xmax)
%% LOGLIKELIHOOD
%
%   Calculates the loglikelihood of measurement scalar
%   varaibility (wm) and produciton scalar variability (wy), given the
%   observed sample (x) and production times (y)

% Use integral
ObsAct = @(m,xmin,xmax,wm)( integral(@(x)(x.*exp(-1./2.*(x-m).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)(exp(-1./2.*(x-m).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true) ./ (1+wy.^2) );
likefunction = @(x,y,m,xmin,xmax,wm,wy)(1./wy./wm./x./ObsAct(m,xmin,xmax,wm).*exp(-1/2.*(y-(ObsAct(m,xmin,xmax,wm))).^2./(wy.*ObsAct(m,xmin,xmax,wm)).^2).*exp(-1/2.*(m-x).^2./(wm.*x).^2));
Likelihoods = integral(@(m)(likefunction(x,y,m,xmin,xmax,wm,wy)),xmin-15*wm*xmin,xmax+15*wm*xmax,'ArrayValued',true);

% Use trapz
% x_vec = min(x):0.1:max(x);
% f_BLS = arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)), m)./arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)./x_vec), m);
% Likelihoods = arrayfun(@(xi,tpi) trapz(1./wy./wm./xi./f_BLS.*exp(-1/2*(tpi-(f_BLS)).^2./(wy.*f_BLS).^2).*exp(-1/2*(m-xi).^2./(wm.*xi).^2)), x, y);

logL = -sum(log(Likelihoods));

function logL = logLikelihoodTRAPZ(wm,wy,m,x,y)
%% LOGLIKELIHOOD
%
%   Calculates the loglikelihood of measurement scalar
%   varaibility (wm) and produciton scalar variability (wy), given the
%   observed sample (x) and production times (y)

% Use trapz
x_vec = min(x):0.1:max(x);
f_BLS = arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)), m)./arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)./x_vec), m);
Likelihoods = arrayfun(@(xi,tpi) trapz(1./wy./wm./xi./f_BLS.*exp(-1/2*(tpi-(f_BLS)).^2./(wy.*f_BLS).^2).*exp(-1/2*(m-xi).^2./(wm.*xi).^2)), x, y);

logL = -sum(log(Likelihoods));

