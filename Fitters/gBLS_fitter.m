function [wm wy b llikelihood] = gBLS_fitter(x,y,varargin)
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
        
        BLSminimizant = @(p)logLikelihoodQUAD(p(1),p(2),p(3),N,x,y,xmin,xmax,dx,M,m);
        
    case 'quad_batch'
        % Use Simpson's quadrature on batches of data
        m = 0:dx:2*xmax;
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
        
        BLSminimizant = @(p)logLikelihoodQUADbatch(p(1),p(2),p(3),N,x,y,xmin,xmax,dx,M,m,batchsize);
        
    case 'quad_nested'        
        % Use Simpson's quadrature on batches of data, performing
        % integrations in nested loops
        m = 0:dx:2*xmax;
        l = length(m);

        BLSminimizant = @(p)logLikelihoodQUADnested(p(1),p(2),p(3),N,x,y,xmin,xmax,dx,m,batchsize);
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
function logL = logLikelihood(wm,wy,g,N,x,y,xmin,xmax)
%% LOGLIKELIHOOD
%
%   Calculates the loglikelihood of measurement scalar
%   varaibility (wm) and produciton scalar variability (wy), given the
%   observed sample (x) and production times (y)

if iscell(N)
    error('Multiple measurement conditions not yet supported for FitType = "integral"')
end

% Use integral
switch N
    case 1
        fBLS = @(m,xmin,xmax,wm)(integral(@(x)(x.*exp(-1./2.*(x-m).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)(exp(-1./2.*(x-m).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true));
        likefunction = @(x,y,m,xmin,xmax,wm,wy,b)(1./wy./wm./x./fBLS(m,xmin,xmax,wm).*exp(-1/2.*(y-(g*fBLS(m,xmin,xmax,wm))).^2./(wy.*g*fBLS(m,xmin,xmax,wm)).^2).*exp(-1/2.*(m-x).^2./(wm.*x).^2));
        Likelihoods = integral(@(m)(likefunction(x,y,m,xmin,xmax,wm,wy,g)),xmin-15*wm*xmin,xmax+15*wm*xmax,'ArrayValued',true);
        
    case 2
        fBLS = @(m1,m2,xmin,xmax,wm)(integral(@(x)(x.*exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)(exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true));
        likefunction = @(x,y,m1,m2,xmin,xmax,wm,wy,g)(1./wy./wm.^2./x.^2./(g*fBLS(m1,m2,xmin,xmax,wm)).*exp(-1/2.*(y-(g*fBLS(m1,m2,xmin,xmax,wm))).^2./(wy.*g*fBLS(m1,m2,xmin,xmax,wm)).^2).*exp(-1/2.*(m1-x).^2./(wm.*x).^2).*exp(-1/2.*(m2-x).^2./(wm.*x).^2));
        for i = 1:length(x)
            Likelihoods(i) = integral2(@(m1,m2)(likefunction(x(i),y(i),m1,m2,xmin,xmax,wm,wy,g)),xmin-15*wm*xmin,xmax+15*wm*xmax,xmin-15*wm*xmin,xmax+15*wm*xmax);
        end
        
end

% Use trapz
% x_vec = min(x):0.1:max(x);
% f_BLS = arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)), m)./arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)./x_vec), m);
% Likelihoods = arrayfun(@(xi,tpi) trapz(1./wy./wm./xi./f_BLS.*exp(-1/2*(tpi-(f_BLS+b)).^2./(wy.*f_BLS).^2).*exp(-1/2*(m-xi).^2./(wm.*xi).^2)), x, y);

logL = -sum(log(Likelihoods));

function logL = logLikelihoodTRAPZ(wm,wy,g,N,m,x,y)
%% LOGLIKELIHOODTRAPZ
%
%   Calculates the loglikelihood of measurement scalar
%   varaibility (wm) and produciton scalar variability (wy), given the
%   observed sample (x) and production times (y)

if iscell(N)
    error('Multiple measurement conditions not yet supported for FitType = "trapz"')
end

% Use trapz
switch N
    case 1
        x_vec = min(x):0.1:max(x);
        f_BLS = arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)), m)./arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)./x_vec), m);
        Likelihoods = arrayfun(@(xi,tpi) trapz(1./wy./wm./xi./(g*f_BLS).*exp(-1/2*(tpi-(g*f_BLS)).^2./(wy.*g*f_BLS).^2).*exp(-1/2*(m-xi).^2./(wm.*xi).^2)), x, y);
        
    case 2
        [M1 M2] = meshgrid(m);
        fBLS = @(m1,m2,xmin,xmax,wm)(integral(@(x)(x.*exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)(exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true));
        integrand = @(m1,m2,xmin,xmax,wm,wp,x,y)( (1./sqrt(2.*pi.*wy.*g*fBLS(m1,m2,xmin,xmax,wm))) .* exp( -(g*fBLS(m1,m2,xmin,xmax,wm)+g-y).^2./(2.*wp.^2.*g*fBLS(m1,m2,xmin,xmax,wm).^2) ) .* (1./sqrt(2.*pi.*wm.^2.*x.^2)).^2 .* exp( -( (m1-x).^2. + (m2-x).^2 )./(2.*wm.^2.*x.^2) ) );
        
        Likelihoods = arrayfun(@(xi,yi) trapz(trapz(integrand(M1,M2,min(x),max(x),wm,wp,xi,yi))), x, y);
        
end

logL = -sum(log(Likelihoods));

function logL = logLikelihoodQUAD(wm,wy,g,N,x,y,xmin,xmax,dx,M,m)
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
%     % Create measument matrix  
%     if xmin-5*wm*xmin < 0
%         m = 0:dx:xmax+5*wm*xmax;
%     else
%         m = xmin-5*wm*xmin:dx:xmax+5*wm*xmax;
%     end
     l = length(m);
%     Mtemp = cell(1,max([N{:}]));
%     [Mtemp{:}] = ndgrid(m);
%     M = zeros(l^max([N{:}]),max([N{:}]));
%     for j = 1:max([N{:}])
%         M(:,j) = [Mtemp{j}(:)];
%     end
    
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
        
        
        method_opts.type = 'quad';
        method_opts.dx = dx;
        fBLS = ScalarBayesEstimators(M(1:l^n,1:n),wm,xmin,xmax,'method',method_opts);
        X = repmat(x{i}',numel(fBLS),1);
        Y = repmat(y{i}',numel(fBLS),1);
        fBLS = repmat(fBLS,1,size(X,2));
        
        p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*g^2*fBLS.^2)) .* exp( -(Y - (g*fBLS)).^2./(2.*wy.^2.*g^2*fBLS.^2) );
        p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*X.^2)).^n .* exp( -sum((repmat(permute(M(1:l^n,1:n),[1 3 2]),[1 size(X,2) 1])-repmat(X,[1 1 n])).^2,3)./(2.*wm.^2.*X.^2) );
        integrand = p_y_take_fBLS.*p_m_take_x;
        
        likelihood = W(1:l^n)'*integrand;
        
        logLi(i) = -sum(log(likelihood));
        
    end
    logL = sum(logLi);

else

%     if xmin-5*wm*xmin < 0
%         m = 0:dx:xmax+5*wm*xmax;
%     else
%         m = xmin-5*wm*xmin:dx:xmax+5*wm*xmax;
%     end
    
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
    
%     Mtemp = cell(1,N);
%     [Mtemp{:}] = ndgrid(m);
%     M = zeros(l^N,N);
%     for i = 1:N
%         M(:,i) = [Mtemp{i}(:)];
%     end
    
    method_opts.type = 'quad';
    method_opts.dx = dx;
    fBLS = ScalarBayesEstimators(M(1:l^N,1:N),wm,xmin,xmax,'method',method_opts);
    X = repmat(x',numel(fBLS),1);
    Y = repmat(y',numel(fBLS),1);
    fBLS = repmat(fBLS,1,size(X,2));
    
    p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*g^2*fBLS.^2)) .* exp( -(Y - (g*fBLS)).^2./(2.*wy.^2.*g^2*fBLS.^2) );
    p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*X.^2)).^N .* exp( -sum((repmat(permute(M(1:l^N,1:N),[1 3 2]),[1 size(X,2) 1])-repmat(X,[1 1 N])).^2,3)./(2.*wm.^2.*X.^2) );
    integrand = p_y_take_fBLS.*p_m_take_x;
    
    likelihood = W(1:l^N)'*integrand;
    
    logL = -sum(log(likelihood));
    
end

function logL = logLikelihoodQUADbatch(wm,wy,g,N,x,y,xmin,xmax,dx,M,m,batchsize)
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
%     % Create measument matrix
%     if xmin-5*wm*xmin < 0
%         m = 0:dx:xmax+5*wm*xmax;
%     else
%         m = xmin-5*wm*xmin:dx:xmax+5*wm*xmax;
%     end
     l = length(m);
%     Mtemp = cell(1,max([N{:}]));
%     [Mtemp{:}] = ndgrid(m);
%     M = zeros(l^max([N{:}]),max([N{:}]));
%     for j = 1:max([N{:}])
%         M(:,j) = [Mtemp{j}(:)];
%     end
    
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
%             if xmin-5*wm*xmin < 0
%                 m = 0:dx:xmax+5*wm*xmax;
%             else
%                 m = xmin-5*wm*xmin:dx:xmax+5*wm*xmax;
%             end
%             
%             if n*length(m)^n > 4000000
%                 error('Surpasing reasonable memory limits; suggest increasing dx or decreasing N')
%             end
            
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

            
            method_opts.type = 'quad';
            method_opts.dx = dx;
            fBLS = ScalarBayesEstimators(M(1:l^n,1:n),wm,xmin,xmax,'method',method_opts);
            X = repmat(xtemp',numel(fBLS),1);
            Y = repmat(ytemp',numel(fBLS),1);
            fBLS = repmat(fBLS,1,size(X,2));
            
            p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*g^2*fBLS.^2)) .* exp( -(Y - (g*fBLS)).^2./(2.*wy.^2.*g^2*fBLS.^2) );
            p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*X.^2)).^n .* exp( -sum((repmat(permute(M(1:l^n,1:n),[1 3 2]),[1 size(X,2) 1])-repmat(X,[1 1 n])).^2,3)./(2.*wm.^2.*X.^2) );
            integrand = p_y_take_fBLS.*p_m_take_x;
            
            likelihood = W(1:l^n)'*integrand;
            
            logLik(i,k) = -sum(log(likelihood));
            
        end
    end
    logL = sum(reshape(logLik,numel(logLik),1));

else
    
%     if xmin-5*wm*xmin < 0
%         m = 0:dx:xmax+5*wm*xmax;
%     else
%         m = xmin-5*wm*xmin:dx:xmax+5*wm*xmax;
%     end
    
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
        
%         % Set up Simpson's nodes
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
        
%         Mtemp = cell(1,N);
%         [Mtemp{:}] = ndgrid(m);
%         M = zeros(l^N,N);
%         for i = 1:N
%             M(:,i) = [Mtemp{i}(:)];
%         end
        
        method_opts.type = 'quad';
        method_opts.dx = dx;
        fBLS = ScalarBayesEstimators(M(1:l^N,1:N),wm,xmin,xmax,'method',method_opts);
        X = repmat(x',numel(fBLS),1);
        Y = repmat(y',numel(fBLS),1);
        fBLS = repmat(fBLS,1,size(X,2));
        
        p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*g^2*fBLS.^2)) .* exp( -(Y - (g*fBLS)).^2./(2.*wy.^2.*g^2*fBLS.^2) );
        p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*X.^2)).^N .* exp( -sum((repmat(permute(M(1:l^N,1:N),[1 3 2]),[1 size(X,2) 1])-repmat(X,[1 1 N])).^2,3)./(2.*wm.^2.*X.^2) );
        integrand = p_y_take_fBLS.*p_m_take_x;
        
        likelihood = W(1:l^N)'*integrand;
        
        logLk(k) = sum(log(likelihood));
    
    end
    logL = -sum(logLk);
    
end

function logL = logLikelihoodQUADnested(wm,wy,g,N,x,y,xmin,xmax,dx,m,batchsize)


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
                    p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*g^2*fBLS.^2)) .* exp( -(ytemp - (g*fBLS)).^2./(2.*wy.^2.*g^2*fBLS.^2) );
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
                        p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*g^2*fBLS.^2)) .* exp( -(ytemp - (g*fBLS)).^2./(2.*wy.^2.*g^2*fBLS.^2) );
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
                            p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*g^2*fBLS.^2)) .* exp( -(ytemp - (g*fBLS)).^2./(2.*wy.^2.*g^2*fBLS.^2) );
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
                                p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*g^2*fBLS.^2)) .* exp( -(ytemp - (g*fBLS)).^2./(2.*wy.^2.*g^2*fBLS.^2) );
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
                p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*g^2*fBLS.^2)) .* exp( -(ytemp - (g*fBLS)).^2./(2.*wy.^2.*g^2*fBLS.^2) );
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
                    p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*g^2*fBLS.^2)) .* exp( -(ytemp - (g*fBLS)).^2./(2.*wy.^2.*g^2*fBLS.^2) );
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
                        p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*g^2*fBLS.^2)) .* exp( -(ytemp - (g*fBLS)).^2./(2.*wy.^2.*g^2*fBLS.^2) );
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
                            p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*g^2*fBLS.^2)) .* exp( -(ytemp - (g*fBLS)).^2./(2.*wy.^2.*g^2*fBLS.^2) );
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