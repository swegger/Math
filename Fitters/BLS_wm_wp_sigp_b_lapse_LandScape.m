function LL = BLS_wm_wp_sigp_b_lapse_LandScape(x,y,parameters,varargin)
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
LapseSupportflg = 0;
CrossValidationflg = 0;
ModelEvidenceflg = 0;
MeasureOptimizationLandscapeflg = 0;
for i = 1:length(varargin)
    if isstr(varargin{i})
        if strcmp(varargin{i},'FitType')
            FitterFlg = 1;
            Fitnum = i;
        elseif strcmp(varargin{i},'N')
            Nflg = 1;
            Nnum = i;
        elseif strcmp(varargin{i},'LapseSupport')
            LapseSupportflg = 1;
            LapseSupportnum = i;
        end
    end
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
            
        case 'quad_batch'
            xmin = varargin{Fitnum+2}(1);
            xmax = varargin{Fitnum+2}(2);
            dx = varargin{Fitnum+2}(3);
            batchsize = varargin{Fitnum+2}(4);
            
        case 'quad_nested'
            xmin = varargin{Fitnum+2}(1);
            xmax = varargin{Fitnum+2}(2);
            dx = varargin{Fitnum+2}(3);
            batchsize = varargin{Fitnum+2}(4);
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

if ~iscell(N)
    if N > 2 & strcmp(FitType,'integral')
        error('N > 2 not yet supported for FitType = "integral"')
    end
end

if LapseSupportflg
    LapseSupport = varargin{LapseSupportnum+1};
else
    if iscell(N)
        LapseSupport = [0 max(vertcat(y{:}))];
    else
        LapseSupport = [0 max(y)];
    end
end


if nargin < 3 
    error('Not enought input arguments!')
end

%% Iterate over different parameter combinations
% Set up likelihood function
switch FitType
    case 'integral'
        % Use integral
        f = @(p)logLikelihood(p(:,1),p(:,2),p(:,3),N,x,y,xmin,xmax);
        
    case 'trapz'
        % Use trapz
        f = @(p)logLikelihoodTRAPZ(p(:,1),p(:,2),p(:,3),N,m,x,y);
        
    case 'quad'
        % Use Simpson's quadrature
        m = 0:dx:2*xmax;
        l = length(m);
        if iscell(N)
            n= max([N{:}]);
            Mtemp = cell(1,n);
            [Mtemp{:}] = ndgrid(m);
            M = zeros(l^n,n);
            for j = 1:n
                M(:,j) = [Mtemp{j}(:)];
            end
        else
            n = N;
            Mtemp = cell(1,N);
            [Mtemp{:}] = ndgrid(m);
            M = zeros(l^N,N);
            for j = 1:N
                M(:,j) = [Mtemp{j}(:)];
            end
        end
        
        f = @(p)logLikelihoodQUAD(p(:,1),p(:,2),p(:,3),p(:,4),p(:,5),N,x,y,xmin,xmax,dx,M,m,LapseSupport(1),LapseSupport(2));
        
    case 'quad_batch'
        % Use Simpson's quadrature on batches of data
        m = 0.05:dx:2*xmax;
        l = length(m);
        if iscell(N)
            n = max([N{:}]);
            Mtemp = cell(1,n);
            [Mtemp{:}] = ndgrid(m);
            M = zeros(l^n,n);
            for j = 1:n
                M(:,j) = [Mtemp{j}(:)];
            end
        else
            n = N;
            Mtemp = cell(1,N);
            [Mtemp{:}] = ndgrid(m);
            M = zeros(l^N,N);
            for j = 1:N
                M(:,j) = [Mtemp{j}(:)];
            end
        end
        
        f = @(p)logLikelihoodQUADbatch(p(:,1),p(:,2),p(:,3),p(:,4),p(:,5),N,x,y,xmin,xmax,dx,M,m,LapseSupport(1),LapseSupport(2),batchsize);
        
    case 'quad_nested'
        % Use Simpson's quadrature on batches of data, performing
        % integrations in nested loops
        m = 0:dx:2*xmax;
        l = length(m);
        
        f = @(p)logLikelihoodQUADnested(p(:,1),p(:,2),p(:,3),N,x,y,xmin,xmax,dx,m,batchsize);
end

for ii = 1:size(parameters,1)
    disp(['Parameter combination ' num2str(ii) ' of ' num2str(size(parameters,1))])
    p = parameters(ii,:);
    LL(ii) = f(p);
    disp(LL(ii))
end    


%% Function to be minimized
function logL = logLikelihood(wm,wy,b,N,x,y,xmin,xmax)
%% LOGLIKELIHOOD
%
%   Calculates the loglikelihood of measurement scalar
%   varaibility (wm) and produciton scalar variability (wy), given the
%   observed sample (x) and production times (y)

error('Lapse model not yet supported for FitType = "integral"')
if iscell(N)
    error('Multiple measurement conditions not yet supported for FitType = "integral"')
end

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

error('Lapse model not yet supported for FitType = "trapz"')
if iscell(N)
    error('Multiple measurement conditions not yet supported for FitType = "trapz"')
end

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

function logL = logLikelihoodQUAD(wm,wy,b,lapse,sig,N,x,y,xmin,xmax,dx,M,m,pmin,pmax)
%% LOGLIKELIHOODQUAD
%
%   Calculates the log likelihood of scalar
%   varaibility (wm) and produciton scalar variability (wy), given the
%   observed sample (x) and production times (y) using Simpson's quadrature
%   method
%
%%



% Determine if different number of Ns are used
if iscell(N)
    l = length(m);
    
    logLi = nan(length(N),length(wm));
    M = repmat(M,[1 1 length(wm)]);
    for i = 1:length(N)
        n = N{i};
        
        if n*length(m)^n > 4000000
            error('Surpasing reasonable memory limits; suggest increasing dx or decreasing N')
        end
        
        % Set up Simpson's nodes
        w = ones(1,l);
        h = (m(end)-m(1))/l;
        w(2:2:l-1) = 4;
        w(3:2:l-1) = 2;
        w = w*h/3;
        
        W = w(:);
        for j = 2:n
            W = W*w;
            W = W(:);
        end
        
        % Lapse Model
        uni = @(pmin,pmax,p)(1/(pmax-pmin));        % Uniform distribution of productions on a lapse trial
        lambda = @(l,p,s)(l);                           % Lapse rate model as a function of sample and production times
        
        % BLS model
        method_opts.type = 'quad';
        method_opts.dx = dx;
        fBLS = nan(size(M(1:l^n,1:n),1),length(wm));
        for ii = 1:length(wm)
            fBLS(:,ii) = ScalarBayesEstimators(M(1:l^n,1:n),wm(ii),xmin,xmax,'method',method_opts);
        end
        X = repmat(x{i}',[size(fBLS,1), 1, length(wm)]);
        Y = repmat(y{i}',[size(fBLS,1), 1, length(wm)]);
        fBLS = repmat(permute(fBLS,[1 3 2]),[1,size(X,2), 1]);
        WM = repmat(permute(wm(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
        WY = repmat(permute(wy(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
        B = repmat(permute(b(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
        SIG = repmat(permute(sig(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
        
        p_y_take_fBLS = (1./sqrt(2.*pi.*(WY.^2.*fBLS.^2 +SIG.^2))) .* exp( -(Y - (fBLS+B)).^2./(2.*(WY.^2.*fBLS.^2 +SIG.^2)) );
        p_m_take_x = (1./sqrt(2.*pi.*WM.^2.*X.^2)).^n .* exp( -squeeze(sum((repmat(permute(M(1:l^n,1:n,:),[1 4 2 3]),[1 size(X,2) 1 1])-repmat(permute(X,[1 2 4 3]),[1 1 n 1])).^2,3))./(2.*WM.^2.*X.^2) );
        integrand = p_y_take_fBLS.*p_m_take_x;
        
        for ii = 1:length(wm)
            likelihood = (1-lambda(lapse(ii),Y,X)).*W(1:l^n)'*integrand(:,:,ii) + lambda(lapse(ii),Y,X).*uni(pmin,pmax,Y(:,:,ii));
            logLi(i,ii) = -sum(log(likelihood),2);
        end
        
    end
    logL = permute(sum(logLi,1),[2 1]);
    
else    
    
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
    for j = 2:N
        W = W*w;
        W = W(:);
    end
    
    % Lapse Model
    uni = @(pmin,pmax,p)(1/(pmax-pmin));        % Uniform distribution of productions on a lapse trial
    lambda = @(l,p,s)(l);                           % Lapse rate model as a function of sample and production times
    
    % BLS model
    method_opts.type = 'quad';
    method_opts.dx = dx;
    fBLS = nan(size(M(1:l^N,1:N),1),length(wm));
    for ii = 1:length(wm)
        fBLS(:,ii) = ScalarBayesEstimators(M(1:l^N,1:N),wm(ii),xmin,xmax,'method',method_opts);
    end
    X = repmat(x',[size(fBLS,1), 1, length(wm)]);
    Y = repmat(y',[size(fBLS,1), 1, length(wm)]);
    M = repmat(M,[1 1 length(wm)]);
    fBLS = repmat(permute(fBLS,[1 3 2]),[1,size(X,2), 1]);
    WM = repmat(permute(wm(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
    WY = repmat(permute(wy(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
    B = repmat(permute(b(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
    SIG = repmat(permute(sig(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
    
    p_y_take_fBLS = (1./sqrt(2.*pi.*(WY.^2.*fBLS.^2 +SIG.^2))) .* exp( -(Y - (fBLS+B)).^2./(2.*(WY.^2.*fBLS.^2 +SIG.^2)) );
    p_m_take_x = (1./sqrt(2.*pi.*WM.^2.*X.^2)).^N .* exp( -squeeze(sum((repmat(permute(M(1:l^N,1:N,:),[1 4 2 3]),[1 size(X,2) 1 1])-repmat(permute(X,[1 2 4 3]),[1 1 N 1])).^2,3))./(2.*WM.^2.*X.^2) );
    integrand = p_y_take_fBLS.*p_m_take_x;
    
    logL = nan(length(wm),1);
    for ii = 1:length(wm)
        likelihood = (1-lambda(lapse(ii),Y,X)).*W(1:l^N)'*integrand(:,:,ii) + lambda(lapse(ii),Y,X).*uni(pmin,pmax,Y(:,:,ii));
        logL(ii,1) = -sum(log(likelihood),2);
    end
    
    
end

function logL = logLikelihoodQUADbatch(wm,wy,b,lapse,sig,N,x,y,xmin,xmax,dx,M,m,pmin,pmax,batchsize)
%% LOGLIKELIHOODQUADbatch
%
%   Calculates the log likelihood of scalar
%   varaibility (wm) and produciton scalar variability (wy), given the
%   observed sample (x) and production times (y) using Simpson's quadrature
%   method
%
%%

% Determine if different number of Ns are used
if iscell(N)

     l = length(m);

for i = 1:length(N)
    logLik{i} = nan(length(wm),ceil(length(x{i})/batchsize));
end
M = repmat(M,[1 1 length(wm)]);
logL = zeros(length(wm),1);
for i = 1:length(N)
    n = N{i};
    
    for k = 1:ceil(length(x{i})/batchsize)
        if length(x{i}) <= (k-1)*batchsize+batchsize
            xtemp = x{i}((k-1)*batchsize+1:end);
            ytemp = y{i}((k-1)*batchsize+1:end);
        else
            xtemp = x{i}((k-1)*batchsize+1:(k-1)*batchsize+batchsize);
            ytemp = y{i}((k-1)*batchsize+1:(k-1)*batchsize+batchsize);
        end
        
        % Set up Simpson's nodes
        w = ones(1,l);
        h = (m(end)-m(1))/l;
        w(2:2:l-1) = 4;
        w(3:2:l-1) = 2;
        w = w*h/3;
        
        W = w(:);
        for j = 2:n
            W = W*w;
            W = W(:);
        end
        
        % Lapse Model
        uni = @(pmin,pmax,p)(1/(pmax-pmin));        % Uniform distribution of productions on a lapse trial
        lambda = @(l,p,s)(l);                           % Lapse rate model as a function of sample and production times
        
        % BLS model
        method_opts.type = 'quad';
        method_opts.dx = dx;
        fBLS = nan(size(M(1:l^n,1:n),1),length(wm));
        for ii = 1:length(wm)
            fBLS(:,ii) = ScalarBayesEstimators(M(1:l^n,1:n),wm(ii),xmin,xmax,'method',method_opts);
        end
        X = repmat(xtemp',[size(fBLS,1), 1, length(wm)]);
        Y = repmat(ytemp',[size(fBLS,1), 1, length(wm)]);
        fBLS = repmat(permute(fBLS,[1 3 2]),[1,size(X,2), 1]);
        WM = repmat(permute(wm(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
        WY = repmat(permute(wy(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
        B = repmat(permute(b(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
        SIG = repmat(permute(sig(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
        
        p_y_take_fBLS = (1./sqrt(2.*pi.*(WY.^2.*fBLS.^2 +SIG.^2))) .* exp( -(Y - (fBLS+B)).^2./(2.*(WY.^2.*fBLS.^2 +SIG.^2)) );
        p_m_take_x = (1./sqrt(2.*pi.*WM.^2.*X.^2)).^n .* exp( -permute(sum((repmat(permute(M(1:l^n,1:n,:),[1 4 2 3]),[1 size(X,2) 1 1])-repmat(permute(X,[1 2 4 3]),[1 1 n 1])).^2,3),[1 2 4 3])./(2.*WM.^2.*X.^2) );
        integrand = p_y_take_fBLS.*p_m_take_x;
        
        for ii = 1:length(wm)
            likelihood = (1-lambda(lapse(ii),Y,X)).*W(1:l^n)'*integrand(:,:,ii) + lambda(lapse(ii),Y,X).*uni(pmin,pmax,Y(:,:,ii));
            logLik{i}(ii,k) = -sum(log(likelihood),2);
        end
    end
    logL = logL + sum(logLik{i},2);
end
    
else

    
    if N*length(m)^N > 4000000
        error('Surpasing reasonable memory limits; suggest increasing dx or decreasing N')
    end
    
    for k = 1:ceil(length(x)/batchsize)
        if length(x) <= (k-1)*batchsize+batchsize
            xtemp = x((k-1)*batchsize+1:end);
            ytemp = y((k-1)*batchsize+1:end);
        else
            xtemp = x((k-1)*batchsize+1:(k-1)*batchsize+batchsize);
            ytemp = y((k-1)*batchsize+1:(k-1)*batchsize+batchsize);
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
        
        
        method_opts.type = 'quad';
        method_opts.dx = dx;
        fBLS = ScalarBayesEstimators(M(1:l^N,1:N),wm,xmin,xmax,'method',method_opts);
        X = repmat(x',numel(fBLS),1);
        Y = repmat(y',numel(fBLS),1);
        fBLS = repmat(fBLS,1,size(X,2));
        
        p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*fBLS.^2)) .* exp( -(Y - (fBLS+b)).^2./(2.*wy.^2.*fBLS.^2) );
        p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*X.^2)).^N .* exp( -sum((repmat(permute(M(1:l^N,1:N),[1 3 2]),[1 size(X,2) 1])-repmat(X,[1 1 N])).^2,3)./(2.*wm.^2.*X.^2) );
        integrand = p_y_take_fBLS.*p_m_take_x;
        
        likelihood = W(1:l^N)'*integrand;
        
        logLk(k) = sum(log(likelihood));
    
    end
    logL = -sum(logLk);
    
end

function logL = logLikelihoodQUADnested(wm,wy,b,N,x,y,xmin,xmax,dx,m,batchsize)


error('Lapse model not yet supported for FitType = "quad_nested"')
% Determine if different number of Ns are used
if iscell(N)
    
    l = length(m);
    
    % For measurement number condition, find likelihood
    for i = 1:length(N)
        n = N{i};
        
        for k = 1:ceil(length(x{i})/batchsize)
            if length(x{i}) <= (k-1)*batchsize+batchsize
                xtemp = x{i}((k-1)*batchsize+1:end);
                ytemp = y{i}((k-1)*batchsize+1:end);
            else
                xtemp = x{i}((k-1)*batchsize+1:(k-1)*batchsize+batchsize);
                ytemp = y{i}((k-1)*batchsize+1:(k-1)*batchsize+batchsize);
            end
            
            
            method_opts.type = 'quad';
            method_opts.dx = dx;
            
            % Used nested loop and serial computation
            logLik(i,k) = 0;
            if n == 1
                % Set up Simpson's nodes
                w = ones(1,l);
                h = (m(end)-m(1))/l;
                w(2:2:l-1) = 4;
                w(3:2:l-1) = 2;
                W = w*h/3;
                integrand = zeros(l,batchsize);
                for ii = 1:l
                    mvec = m(ii);
                    fBLS = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                    p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*fBLS.^2)) .* exp( -(ytemp - (fBLS+b)).^2./(2.*wy.^2.*fBLS.^2) );
                    p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                    integrand(ii,:) = p_y_take_fBLS.*p_m_take_x;
                    ind = ii;
                end
                likelihood = W*integrand; %sum(repmat(w',1,size(integrand,2)).*integrand);
                logLik(i,k) = logLik(i,k) - sum(log(likelihood));
                
            elseif n == 2
                % Set up Simpson's nodes
                w = ones(1,l);
                h = (m(end)-m(1))/l;
                w(2:2:l-1) = 4;
                w(3:2:l-1) = 2;
                w = w*h/3;
                
                W = w(:);
                W = W*w;
                W = W(:);
                integrand = zeros(l^n,batchsize);
                for ii = 1:l
                    for jj = 1:l
                        mvec = [m(ii) m(jj)];
                        fBLS = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                        p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*fBLS.^2)) .* exp( -(ytemp - (fBLS+b)).^2./(2.*wy.^2.*fBLS.^2) );
                        p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                        ind = l*(jj-1) + ii;
                        integrand(ind,:) = p_y_take_fBLS.*p_m_take_x;
                    end
                end
                likelihood = W'*integrand;
                logLik(i,k) = logLik(i,k) - sum(log(likelihood));
                
            elseif n == 3
                % Set up Simpson's nodes
                w = ones(1,l);
                h = (m(end)-m(1))/l;
                w(2:2:l-1) = 4;
                w(3:2:l-1) = 2;
                w = w*h/3;
                
                W = w(:);
                for j = 2:n
                    W = W*w;
                    W = W(:);
                end
                integrand = zeros(l^n,batchsize);
                for ii = 1:l
                    for jj = 1:l
                        for kk = 1:l
                            mvec = [m(ii) m(jj) m(kk)];
                            fBLS = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                            p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*fBLS.^2)) .* exp( -(ytemp - (fBLS+b)).^2./(2.*wy.^2.*fBLS.^2) );
                            p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                            ind = l^2*(kk-1) + l*(jj-1) + ii;
                            integrand(ind,:) = p_y_take_fBLS.*p_m_take_x;
                        end
                    end
                end
                likelihood = W'*integrand;
                logLik(i,k) = logLik(i,k) - sum(log(likelihood));
                
            elseif n == 4
                % Set up Simpson's nodes
                w = ones(1,l);
                h = (m(end)-m(1))/l;
                w(2:2:l-1) = 4;
                w(3:2:l-1) = 2;
                w = w*h/3;
                
                W = w(:);
                for j = 2:n
                    W = W*w;
                    W = W(:);
                end
                
                integrand = zeros(l^n,batchsize);
                for ii = 1:l
                    for jj = 1:l
                        for kk = 1:l
                            for ll = 1:l
                                mvec = [m(ii) m(jj) m(kk) m(ll)];
                                fBLS = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                                p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*fBLS.^2)) .* exp( -(ytemp - (fBLS+b)).^2./(2.*wy.^2.*fBLS.^2) );
                                p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                                ind = l^3*(ll-1) + l^2*(kk-1) + l*(jj-1) + ii;
                                integrand(ind,:) = p_y_take_fBLS.*p_m_take_x;
                            end
                        end
                    end
                end
                likelihood = W'*integrand;
                logLik(i,k) = logLik(i,k) - sum(log(likelihood));
            else
                error('Only 1-4 measurements supported')
            end
            
            
        end
    end
    logL = sum(reshape(logLik,numel(logLik),1));
    
else
    
    for k = 1:ceil(length(x)/batchsize)
        if length(x) <= (k-1)*batchsize+batchsize
            xtemp = x((k-1)*batchsize+1:end);
            ytemp = y((k-1)*batchsize+1:end);
        else
            xtemp = x((k-1)*batchsize+1:(k-1)*batchsize+batchsize);
            ytemp = y((k-1)*batchsize+1:(k-1)*batchsize+batchsize);
        end
        
        
        method_opts.type = 'quad';
        method_opts.dx = dx;
        
        % Used nested loop and serial computation
        logLk(k) = 0;
        if n == 1
            % Set up Simpson's nodes
            w = ones(1,l);
            h = (m(end)-m(1))/l;
            w(2:2:l-1) = 4;
            w(3:2:l-1) = 2;
            W = w*h/3;
            integrand = zeros(l,batchsize);
            for ii = 1:l
                mvec = m(ii);
                fBLS = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*fBLS.^2)) .* exp( -(ytemp - (fBLS+b)).^2./(2.*wy.^2.*fBLS.^2) );
                p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                integrand(ii,:) = p_y_take_fBLS.*p_m_take_x;
                ind = ii;
            end
            likelihood = W*integrand; %sum(repmat(w',1,size(integrand,2)).*integrand);
            logLik(k) = logLik(k) - sum(log(likelihood));
            
        elseif n == 2
            % Set up Simpson's nodes
            w = ones(1,l);
            h = (m(end)-m(1))/l;
            w(2:2:l-1) = 4;
            w(3:2:l-1) = 2;
            w = w*h/3;
            
            W = w(:);
            W = W*w;
            W = W(:);
            integrand = zeros(l^n,batchsize);
            for ii = 1:l
                for jj = 1:l
                    mvec = [m(ii) m(jj)];
                    fBLS = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                    p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*fBLS.^2)) .* exp( -(ytemp - (fBLS+b)).^2./(2.*wy.^2.*fBLS.^2) );
                    p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                    ind = l*(jj-1) + ii;
                    integrand(ind,:) = p_y_take_fBLS.*p_m_take_x;
                end
            end
            likelihood = W'*integrand;
            logLik(k) = logLik(k) - sum(log(likelihood));
            
        elseif n == 3
            % Set up Simpson's nodes
            w = ones(1,l);
            h = (m(end)-m(1))/l;
            w(2:2:l-1) = 4;
            w(3:2:l-1) = 2;
            w = w*h/3;
            
            W = w(:);
            for j = 2:n
                W = W*w;
                W = W(:);
            end
            integrand = zeros(l^n,batchsize);
            for ii = 1:l
                for jj = 1:l
                    for kk = 1:l
                        mvec = [m(ii) m(jj) m(kk)];
                        fBLS = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                        p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*fBLS.^2)) .* exp( -(ytemp - (fBLS+b)).^2./(2.*wy.^2.*fBLS.^2) );
                        p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                        ind = l^2*(kk-1) + l*(jj-1) + ii;
                        integrand(ind,:) = p_y_take_fBLS.*p_m_take_x;
                    end
                end
            end
            likelihood = W'*integrand;
            logLik(k) = logLik(k) - sum(log(likelihood));
            
        elseif n == 4
            % Set up Simpson's nodes
            w = ones(1,l);
            h = (m(end)-m(1))/l;
            w(2:2:l-1) = 4;
            w(3:2:l-1) = 2;
            w = w*h/3;
            
            W = w(:);
            for j = 2:n
                W = W*w;
                W = W(:);
            end
            
            integrand = zeros(l^n,batchsize);
            for ii = 1:l
                for jj = 1:l
                    for kk = 1:l
                        for ll = 1:l
                            mvec = [m(ii) m(jj) m(kk) m(ll)];
                            fBLS = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                            p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*fBLS.^2)) .* exp( -(ytemp - (fBLS+b)).^2./(2.*wy.^2.*fBLS.^2) );
                            p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                            ind = l^3*(ll-1) + l^2*(kk-1) + l*(jj-1) + ii;
                            integrand(ind,:) = p_y_take_fBLS.*p_m_take_x;
                        end
                    end
                end
            end
            likelihood = W'*integrand;
            logLik(k) = logLik(k) - sum(log(likelihood));
        else
            error('Only 1-4 measurements supported')
        end
        
        
    end
    logL = -sum(logLk);
    
end


%% ModelEvidence functions
    function out = modelEvidenceFunction(ll,lfun)
        out = exp(-lfun + ll);
        out(isnan(out)) = 0;
        out = out(:);