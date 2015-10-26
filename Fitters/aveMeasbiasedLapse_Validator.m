function [llikelihood, lmodelEvidence] = aveMeasbiasedLapse_Validator(x,y,varargin)
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
ModelEvidenceflg = 0;
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
        elseif strcmp(varargin{i},'LapseSupport')
            LapseSupportflg = 1;
            LapseSupportnum = i;
        elseif strcmp(varargin{i},'ModelEvidence')
            ModelEvidenceflg = 1;
            ModelEvidencenum = i;
        end
    end
end

if ICflg
    IC = varargin{ICnum+1};
else
    IC = [0.1 0.06 0 0.05];
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

if ModelEvidenceflg
    ModelEvidence = varargin{ModelEvidencenum+1};
    if ~isfield(ModelEvidence,'OpenMind')
        ModelEvidence.OpenMind = 0;
    end
else
    ModelEvidence.method = 'none';
end

if nargin < 3 
    error('Not enought input arguments!')
end

%% Test model predictions on data
for ii = 1:length(x)
    switch FitType
        case 'integral'
            % Use integral
            validant = @(p)logLikelihood(p(:,1),p(:,2),p(:,3),N,x,y,xmin,xmax);
            
        case 'trapz'
            % Use trapz
            validant = @(p)logLikelihoodTRAPZ(p(:,1),p(:,2),p(:,3),N,m,x,y);
            
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
            
            validant = @(p)logLikelihoodQUAD(p(:,1),p(:,2),p(:,3),p(:,4),N,x,y,xmin,xmax,dx,M,m,LapseSupport(1),LapseSupport(2));
            
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
            
            validant = @(p)logLikelihoodQUADbatch(p(:,1),p(:,2),p(:,3),p(:,4),N,x,y,xmin,xmax,dx,M,m,LapseSupport(1),LapseSupport(2),batchsize);
            
        case 'quad_nested'
            % Use Simpson's quadrature on batches of data, performing
            % integrations in nested loops
            m = 0:dx:2*xmax;
            l = length(m);
            
            validant = @(p)logLikelihoodQUADnested(p(:,1),p(:,2),p(:,3),N,x,y,xmin,xmax,dx,m,batchsize);
    end
    lparams = [wm, wp, b, lapse];
    llikelihood(ii) = validant(lparams);
            
    
    switch ModelEvidence.method
        case 'Integration'
            % Numerical integration of the likelihood function across
            % different parameter values.
            functionHandle = @(P)(modelEvidenceFunction(llikelihood(ii),validant(P)));
            modelEvidence = ndintegrate(functionHandle,ModelEvidence.paramLimits,'method',ModelEvidence.integrationMethod,'options',ModelEvidence.integrationOptions,'OpenMind',ModelEvidence.OpenMind);
            lmodelEvidence(ii) = log(modelEvidence) - llikelihood(ii);      % log model evidence
            
        case 'UniformApproximation'
            error('UniformApproximation method for model evidence not yet supported!')
            
        case 'GaussianApproximation'
            % Use Laplace's method to estimate the integral over the
            % parameters via a Gaussian with mu at the ML fit and
            % covariance from the Hessian at the ML fit.
            error('Gaussian approximation method for model evidence not yet supported!')
            %lmodelEvidence(ii) = -llikelihood(ii) + log(1/prod(diff(ModelEvidence.paramLimits,1,2))) + length(lparams)/2*log(2*pi) - log(det(-hessian))/2;  % log model evidence <- this hessian isn't right... its the hessian of ln(p(D|p,M))
            
        case {'None','none'}
            lmodelEvidence(ii) = NaN;
            
        otherwise
            error(['Model evidence calculation method ' ModelEvidence.method ' not recognized!'])
    end

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
        f = @(m,xmin,xmax,wm)(integral(@(x)(x.*exp(-1./2.*(x-m).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)(exp(-1./2.*(x-m).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true));
        likefunction = @(x,y,m,xmin,xmax,wm,wy,b)(1./wy./wm./x./f(m,xmin,xmax,wm).*exp(-1/2.*(y-(f(m,xmin,xmax,wm)+b)).^2./(wy.*f(m,xmin,xmax,wm)).^2).*exp(-1/2.*(m-x).^2./(wm.*x).^2));
        Likelihoods = integral(@(m)(likefunction(x,y,m,xmin,xmax,wm,wy,b)),xmin-15*wm*xmin,xmax+15*wm*xmax,'ArrayValued',true);
        
    case 2
        f = @(m1,m2,xmin,xmax,wm)(integral(@(x)(x.*exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)(exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true));
        likefunction = @(x,y,m1,m2,xmin,xmax,wm,wy,b)(1./wy./wm.^2./x.^2./f(m1,m2,xmin,xmax,wm).*exp(-1/2.*(y-(f(m1,m2,xmin,xmax,wm)+b)).^2./(wy.*f(m1,m2,xmin,xmax,wm)).^2).*exp(-1/2.*(m1-x).^2./(wm.*x).^2).*exp(-1/2.*(m2-x).^2./(wm.*x).^2));
        for i = 1:length(x)
            Likelihoods(i) = integral2(@(m1,m2)(likefunction(x(i),y(i),m1,m2,xmin,xmax,wm,wy,b)),xmin-15*wm*xmin,xmax+15*wm*xmax,xmin-15*wm*xmin,xmax+15*wm*xmax);
        end
        
end

% Use trapz
% x_vec = min(x):0.1:max(x);
% f_ = arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)), m)./arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)./x_vec), m);
% Likelihoods = arrayfun(@(xi,tpi) trapz(1./wy./wm./xi./f_.*exp(-1/2*(tpi-(f_+b)).^2./(wy.*f_).^2).*exp(-1/2*(m-xi).^2./(wm.*xi).^2)), x, y);

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
        f_ = arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)), m)./arrayfun(@(mi) trapz(exp(-1./2.*(x_vec-mi).^2./wm.^2./x_vec.^2)./x_vec), m);
        Likelihoods = arrayfun(@(xi,tpi) trapz(1./wy./wm./xi./f_.*exp(-1/2*(tpi-(f_+b)).^2./(wy.*f_).^2).*exp(-1/2*(m-xi).^2./(wm.*xi).^2)), x, y);
        
    case 2
        [M1 M2] = meshgrid(m);
        f = @(m1,m2,xmin,xmax,wm)(integral(@(x)(x.*exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)(exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true));
        integrand = @(m1,m2,xmin,xmax,wm,wp,x,y)( (1./sqrt(2.*pi.*wy.*f(m1,m2,xmin,xmax,wm))) .* exp( -(f(m1,m2,xmin,xmax,wm)+b-y).^2./(2.*wp.^2.*f(m1,m2,xmin,xmax,wm).^2) ) .* (1./sqrt(2.*pi.*wm.^2.*x.^2)).^2 .* exp( -( (m1-x).^2. + (m2-x).^2 )./(2.*wm.^2.*x.^2) ) );
        
        Likelihoods = arrayfun(@(xi,yi) trapz(trapz(integrand(M1,M2,min(x),max(x),wm,wp,xi,yi))), x, y);
        
end

logL = -sum(log(Likelihoods));

function logL = logLikelihoodQUAD(wm,wy,b,lapse,N,x,y,xmin,xmax,dx,M,m,pmin,pmax)
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
        
        % Lapse Model
        uni = @(pmin,pmax,p)(1/(pmax-pmin));        % Uniform distribution of productions on a lapse trial
        lambda = @(l,p,s)(l);                           % Lapse rate model as a function of sample and production times
        
        % Averaging model
        method_opts.type = 'quad';
        method_opts.dx = dx;
        estimator.type = 'weightedMean';
        estimator.weights = ones(1,n)/n;
        f = ScalarBayesEstimators(M(1:l^n,1:n),wm,xmin,xmax,'method',method_opts,'estimator',estimator);
        X = repmat(x{i}',numel(f),1);
        Y = repmat(y{i}',numel(f),1);
        f = repmat(f,1,size(X,2));
        
        p_y_take_f = (1./sqrt(2.*pi.*wy.^2.*f.^2)) .* exp( -(Y - (f+b)).^2./(2.*wy.^2.*f.^2) );
        p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*X.^2)).^n .* exp( -sum((repmat(permute(M(1:l^n,1:n),[1 3 2]),[1 size(X,2) 1])-repmat(X,[1 1 n])).^2,3)./(2.*wm.^2.*X.^2) );
%        p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*X.^2)) .* exp( -(repmat(permute(sum(M(1:l^n,1:n).*repmat(estimator.weights,l^n,1),2),[1 3 2]),[1 size(X,2) 1])-X).^2./(2.*wm.^2.*X.^2) );
        integrand = p_y_take_f.*p_m_take_x;
        
        likelihood = (1-lambda(lapse,Y,X)).*W(1:l^n)'*integrand + lambda(lapse,Y,X).*uni(pmin,pmax,Y);
        
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
    

    % Lapse Model
    uni = @(pmin,pmax,p)(1/(pmax-pmin));        % Uniform distribution of productions on a lapse trial
    lambda = @(l,p,s)(l);                           % Lapse rate model as a function of sample and production times
    
    % Avering model
    method_opts.type = 'quad';
    method_opts.dx = dx;
    f = ScalarBayesEstimators(M(1:l^N,1:N),wm,xmin,xmax,'method',method_opts);
    X = repmat(x',numel(f),1);
    Y = repmat(y',numel(f),1);
    f = repmat(f,1,size(X,2));
    
    p_y_take_f = (1./sqrt(2.*pi.*wy.^2.*f.^2)) .* exp( -(Y - (f+b)).^2./(2.*wy.^2.*f.^2) );
    p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*X.^2)).^N .* exp( -sum((repmat(permute(M(1:l^N,1:N),[1 3 2]),[1 size(X,2) 1])-repmat(X,[1 1 N])).^2,3)./(2.*wm.^2.*X.^2) );
    integrand = p_y_take_f.*p_m_take_x;
    
    likelihood = (1-lambda(lapse,Y,X)).*W(1:l^N)'*integrand + lambda(lapse,Y,X)*uni(pmin,pmax,Y);
    
    logL = -sum(log(likelihood));
    
end

function logL = logLikelihoodQUADbatch(wm,wy,b,lapse,N,x,y,xmin,xmax,dx,M,m,pmin,pmax,batchsize)
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
    
%     % For measurement number condition, find likelihood
%     for i = 1:length(N)
%         n = N{i};
%         
%         for k = 1:ceil(length(x{i})/batchsize)
%             if length(x{i}) <= (k-1)*batchsize+batchsize
%                 xtemp = x{i}((k-1)*batchsize+1:end);
%                 ytemp = y{i}((k-1)*batchsize+1:end);
%             else
%                 xtemp = x{i}((k-1)*batchsize+1:(k-1)*batchsize+batchsize);
%                 ytemp = y{i}((k-1)*batchsize+1:(k-1)*batchsize+batchsize);
%             end
% %             if xmin-5*wm*xmin < 0
% %                 m = 0:dx:xmax+5*wm*xmax;
% %             else
% %                 m = xmin-5*wm*xmin:dx:xmax+5*wm*xmax;
% %             end
% %             
% %             if n*length(m)^n > 4000000
% %                 error('Surpasing reasonable memory limits; suggest increasing dx or decreasing N')
% %             end
%             
%             % Set up Simpson's nodes
%             w = ones(1,l);
%             h = (m(end)-m(1))/l;
%             w(2:2:l-1) = 4;
%             w(3:2:l-1) = 2;
%             w = w*h/3;
%             
%             W = w(:);
%             for j = 2:n
%                 W = W*w;
%                 W = W(:);
%             end
% 
%             
%             method_opts.type = 'quad';
%             method_opts.dx = dx;
%             f = ScalarBayesEstimators(M(1:l^n,1:n),wm,xmin,xmax,'method',method_opts);
%             X = repmat(xtemp',numel(f),1);
%             Y = repmat(ytemp',numel(f),1);
%             f = repmat(f,1,size(X,2));
%             
%             p_y_take_f = (1./sqrt(2.*pi.*wy.^2.*f.^2)) .* exp( -(Y - (f+b)).^2./(2.*wy.^2.*f.^2) );
%             p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*X.^2)).^n .* exp( -sum((repmat(permute(M(1:l^n,1:n),[1 3 2]),[1 size(X,2) 1])-repmat(X,[1 1 n])).^2,3)./(2.*wm.^2.*X.^2) );
%             integrand = p_y_take_f.*p_m_take_x;
%             
%             likelihood = W(1:l^n)'*integrand;
%             
%             logLik(i,k) = -sum(log(likelihood));
%             
%         end
%     end
%     logL = sum(reshape(logLik,numel(logLik),1));

    
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
            estimator.type = 'weightedMean';
            estimator.weights = ones(1,n)/n;
            fBLS = nan(size(M(1:l^n,1:n),1),length(wm));
            for ii = 1:length(wm)
                fBLS(:,ii) = ScalarBayesEstimators(M(1:l^n,1:n),wm(ii),xmin,xmax,'method',method_opts,'estimator',estimator);
            end
            X = repmat(xtemp',[size(fBLS,1), 1, length(wm)]);
            Y = repmat(ytemp',[size(fBLS,1), 1, length(wm)]);
            fBLS = repmat(permute(fBLS,[1 3 2]),[1,size(X,2), 1]);
            WM = repmat(permute(wm(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
            WY = repmat(permute(wy(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
            B = repmat(permute(b(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
            %LAPSE = repmat(permute(lapse(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
            
            p_y_take_fBLS = (1./sqrt(2.*pi.*WY.^2.*fBLS.^2)) .* exp( -(Y - (fBLS+B)).^2./(2.*WY.^2.*fBLS.^2) );
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
        f = ScalarBayesEstimators(M(1:l^N,1:N),wm,xmin,xmax,'method',method_opts);
        X = repmat(x',numel(f),1);
        Y = repmat(y',numel(f),1);
        f = repmat(f,1,size(X,2));
        
        p_y_take_f = (1./sqrt(2.*pi.*wy.^2.*f.^2)) .* exp( -(Y - (f+b)).^2./(2.*wy.^2.*f.^2) );
        p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*X.^2)).^N .* exp( -sum((repmat(permute(M(1:l^N,1:N),[1 3 2]),[1 size(X,2) 1])-repmat(X,[1 1 N])).^2,3)./(2.*wm.^2.*X.^2) );
        integrand = p_y_take_f.*p_m_take_x;
        
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
                    f = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                    p_y_take_f = (1./sqrt(2.*pi.*wy.^2.*f.^2)) .* exp( -(ytemp - (f+b)).^2./(2.*wy.^2.*f.^2) );
                    p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                    integrand(ii,:) = p_y_take_f.*p_m_take_x;
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
                        f = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                        p_y_take_f = (1./sqrt(2.*pi.*wy.^2.*f.^2)) .* exp( -(ytemp - (f+b)).^2./(2.*wy.^2.*f.^2) );
                        p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                        ind = l*(jj-1) + ii;
                        integrand(ind,:) = p_y_take_f.*p_m_take_x;
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
                            f = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                            p_y_take_f = (1./sqrt(2.*pi.*wy.^2.*f.^2)) .* exp( -(ytemp - (f+b)).^2./(2.*wy.^2.*f.^2) );
                            p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                            ind = l^2*(kk-1) + l*(jj-1) + ii;
                            integrand(ind,:) = p_y_take_f.*p_m_take_x;
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
                                f = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                                p_y_take_f = (1./sqrt(2.*pi.*wy.^2.*f.^2)) .* exp( -(ytemp - (f+b)).^2./(2.*wy.^2.*f.^2) );
                                p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                                ind = l^3*(ll-1) + l^2*(kk-1) + l*(jj-1) + ii;
                                integrand(ind,:) = p_y_take_f.*p_m_take_x;
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
                f = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                p_y_take_f = (1./sqrt(2.*pi.*wy.^2.*f.^2)) .* exp( -(ytemp - (f+b)).^2./(2.*wy.^2.*f.^2) );
                p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                integrand(ii,:) = p_y_take_f.*p_m_take_x;
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
                    f = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                    p_y_take_f = (1./sqrt(2.*pi.*wy.^2.*f.^2)) .* exp( -(ytemp - (f+b)).^2./(2.*wy.^2.*f.^2) );
                    p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                    ind = l*(jj-1) + ii;
                    integrand(ind,:) = p_y_take_f.*p_m_take_x;
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
                        f = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                        p_y_take_f = (1./sqrt(2.*pi.*wy.^2.*f.^2)) .* exp( -(ytemp - (f+b)).^2./(2.*wy.^2.*f.^2) );
                        p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                        ind = l^2*(kk-1) + l*(jj-1) + ii;
                        integrand(ind,:) = p_y_take_f.*p_m_take_x;
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
                            f = ScalarBayesEstimators(mvec,wm,xmin,xmax,'method',method_opts);
                            p_y_take_f = (1./sqrt(2.*pi.*wy.^2.*f.^2)) .* exp( -(ytemp - (f+b)).^2./(2.*wy.^2.*f.^2) );
                            p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*xtemp.^2)).^n .* exp( -sum((repmat(mvec,length(xtemp),1)-repmat(xtemp(:),1,n)).^2,2)./(2.*wm.^2.*xtemp.^2) );
                            ind = l^3*(ll-1) + l^2*(kk-1) + l*(jj-1) + ii;
                            integrand(ind,:) = p_y_take_f.*p_m_take_x;
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