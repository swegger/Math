function [wm wy b llikelihood] = BLSbiased_fitter(x,y,varargin)
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
Nflg = 0;
for i = 1:length(varargin)
    if isstr(varargin{i})
        if strcmp(varargin{i},'InitCond')
            ICflg = 1;
            ICnum = i;
        elseif strcmp(varargin{i},'FitType')
            FitterFlg = 1;
            Fitnum = i;
        elseif strcmp(varargin{i},'N')
            Nflg = 1;
            Nnum = i;
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
            
        case 'quad'
            xmin = varargin{Fitnum+2}(1);
            xmax = varargin{Fitnum+2}(2);
            dx = varargin{Fitnum+2}(3);
    end
        
else
    FitType = 'integral';
    xmin = min(x);
    xmax = max(x);
end

if Nflg
    N = varargin{Nnum+1};
else
    N = 1;
end

if N > 2
    error('N > 2 not yet supported')
end


if nargin < 3
    error('Not enought input arguments!')
end


% Fit according to minimizer
OPTIONS = optimset('Display','iter');

switch FitType
    case 'integral'
        % Use integral
        BLSminimizant = @(p)logLikelihood(p(1),p(2),p(3),N,x,y,xmin,xmax);
        
    case 'trapz'
        % Use trapz
        BLSminimizant = @(p)logLikelihoodTRAPZ(p(1),p(2),p(3),N,m,x,y);
        
    case 'quad'
        % Use Simpson's quadrature
        BLSminimizant = @(p)logLikelihoodQUAD(p(1),p(2),p(3),N,x,y,xmin,xmax,dx);
        
end

minimizer = 'fminsearch(BLSminimizant, [wM_ini wP_ini b_ini], OPTIONS);';
wM_ini = IC(1);
wP_ini = IC(2);
b_ini = IC(3);
[lparams llikelihood] = eval(minimizer);

wm = lparams(1);
wy = lparams(2);
b = lparams(3);
llikelihood = llikelihood;
    

%% Function to be minimized
function logL = logLikelihood(wm,wy,b,N,x,y,xmin,xmax)
%% LOGLIKELIHOOD
%
%   Calculates the loglikelihood of measurement scalar
%   varaibility (wm) and produciton scalar variability (wy), given the
%   observed sample (x) and production times (y)

% Use integral
switch N
    case 1
        fBLS = @(m,xmin,xmax,wm)(integral(@(x)(x.*exp(-1./2.*(x-m).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)(exp(-1./2.*(x-m).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true));
        likefunction = @(x,y,m,xmin,xmax,wm,wy,b)(1./wy./wm./x./fBLS(m,xmin,xmax,wm).*exp(-1/2.*(y-(fBLS(m,xmin,xmax,wm)+b)).^2./(wy.*fBLS(m,xmin,xmax,wm)).^2).*exp(-1/2.*(m-x).^2./(wm.*x).^2));
        Likelihoods = integral(@(m)(likefunction(x,y,m,xmin,xmax,wm,wy,b)),xmin-15*wm*xmin,xmax+15*wm*xmax,'ArrayValued',true);
        
    case 2
        fBLS = @(m1,m2,xmin,xmax,wm)(integral(@(x)(x.*exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)(exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true));
        likefunction = @(x,y,m1,m2,xmin,xmax,wm,wy,b)(1./wy./wm.^2./x.^2./fBLS(m1,m2,xmin,xmax,wm).*exp(-1/2.*(y-(fBLS(m1,m2,xmin,xmax,wm)+b)).^2./(wy.*fBLS(m1,m2,xmin,xmax,wm)).^2).*exp(-1/2.*(m1-x).^2./(wm.*x).^2).*exp(-1/2.*(m2-x).^2./(wm.*x).^2));
        for i = 1:length(x)
            Likelihoods(i) = integral2(@(m1,m2)(likefunction(x(i),y(i),m1,m2,xmin,xmax,wm,wy,b)),xmin-15*wm*xmin,xmax+15*wm*xmax,xmin-15*wm*xmin,xmax+15*wm*xmax);
        end
        
end

% Use trapz
% x_vec = min(x):0.1:max(x);
% f_BLS = arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)), m)./arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)./x_vec), m);
% Likelihoods = arrayfun(@(xi,tpi) trapz(1./wy./wm./xi./f_BLS.*exp(-1/2*(tpi-(f_BLS+b)).^2./(wy.*f_BLS).^2).*exp(-1/2*(m-xi).^2./(wm.*xi).^2)), x, y);

logL = -sum(log(Likelihoods));

function logL = logLikelihoodTRAPZ(wm,wy,b,N,m,x,y)
%% LOGLIKELIHOODTRAPZ
%
%   Calculates the loglikelihood of measurement scalar
%   varaibility (wm) and produciton scalar variability (wy), given the
%   observed sample (x) and production times (y)

% Use trapz
switch N
    case 1
        x_vec = min(x):0.1:max(x);
        f_BLS = arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)), m)./arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)./x_vec), m);
        Likelihoods = arrayfun(@(xi,tpi) trapz(1./wy./wm./xi./f_BLS.*exp(-1/2*(tpi-(f_BLS+b)).^2./(wy.*f_BLS).^2).*exp(-1/2*(m-xi).^2./(wm.*xi).^2)), x, y);
        
    case 2
        [M1 M2] = meshgrid(m);
        fBLS = @(m1,m2,xmin,xmax,wm)(integral(@(x)(x.*exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)(exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true));
        integrand = @(m1,m2,xmin,xmax,wm,wp,x,y)( (1./sqrt(2.*pi.*wy.*fBLS(m1,m2,xmin,xmax,wm))) .* exp( -(fBLS(m1,m2,xmin,xmax,wm)+b-y).^2./(2.*wp.^2.*fBLS(m1,m2,xmin,xmax,wm).^2) ) .* (1./sqrt(2.*pi.*wm.^2.*x.^2)).^2 .* exp( -( (m1-x).^2. + (m2-x).^2 )./(2.*wm.^2.*x.^2) ) );
        
        Likelihoods = arrayfun(@(xi,yi) trapz(trapz(integrand(M1,M2,min(x),max(x),wm,wp,xi,yi))), x, y);
        
end

logL = -sum(log(Likelihoods));

function logL = logLikelihoodQUAD(wm,wy,b,N,x,y,xmin,xmax,dx)
%% LOGLIKELIHOODQUAD
%
%   Calculates the log likelihood of scalar
%   varaibility (wm) and produciton scalar variability (wy), given the
%   observed sample (x) and production times (y) using Simpson's quadrature
%   method
%
%%

% use quadradic interpolation (Simpson's method)
if xmin-5*wm*xmin < 0
    m = 0:dx:xmax+5*wm*xmax;
else
    m = xmin-5*wm*xmin:dx:xmax+5*wm*xmax;
end

if N*length(m)^N > 4000000
    error('Surpasing reasonable memory limits; suggest increasing dx or decreasing N')
end

% Set up Simpson's nodes
l = length(m);
w = ones(1,l);
h = (m(end)-m(1))/l;
w(2:2:l-1) = 4;
w(3:2:l-1) = 2;
w = w*h/3;

W = w(:);
for i = 2:N
    W = W*w;
    W = W(:);
end

Mtemp = cell(1,N);
[Mtemp{:}] = ndgrid(m);
for i = 1:N
    M(:,i) = [Mtemp{i}(:)];
end

method_opts.type = 'quad';
method_opts.dx = dx;
fBLS = ScalarBayesEstimators(M,wm,xmin,xmax,'method',method_opts);
X = repmat(x',numel(fBLS),1);
Y = repmat(y',numel(fBLS),1);
fBLS = repmat(fBLS,1,size(X,2));

p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*fBLS.^2)) .* exp( -(Y - (fBLS+b)).^2./(2.*wy.^2.*fBLS.^2) );
p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*X.^2)).^N .* exp( -sum((repmat(permute(M,[1 3 2]),[1 size(X,2) 1])-repmat(X,[1 1 N])).^2,3)./(2.*wm.^2.*X.^2) );
integrand = p_y_take_fBLS.*p_m_take_x;

likelihood = W'*integrand;

logL = -sum(log(likelihood));