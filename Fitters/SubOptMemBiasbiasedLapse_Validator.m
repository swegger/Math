function [llikelihood, lmodelEvidence] = SubOptMemBiasbiasedLapse_Validator(x,y,wm,wy,wm_drift,w_int,b,lapse,varargin)
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
    IC = [0.1 0.06 0.1 0.1/sqrt(2) 0 0.05];
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
            error('Not yet supported!')
            
        case 'trapz'
            error('Not yet supported!')
            
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
            
            validant = @(p)logLikelihoodQUAD(p(:,1),p(:,2),p(:,3),p(:,4),p(:,5),p(:,6),N,x,y,xmin,xmax,dx,M,m,LapseSupport(1),LapseSupport(2));
            
        case 'quad_batch'
            error('Not yet supported!')
            
        case 'quad_nested'
            error('Not yet supported!')
    end
    lparams = [wm, wy, wm_drift, w_int, b, lapse];
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

function logL = logLikelihoodQUAD(wm,wy,b,wm_drift,w_int,lapse,N,x,y,xmin,xmax,dx,M,m,pmin,pmax)
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
        estimator.type = 'SubOptMemBias';
        estimator.wm_drift = wm_drift;
        estimator.w_int = w_int;
        fBLS = nan(size(M(1:l^n,1:n),1),length(wm));
        for ii = 1:length(wm)
            fBLS(:,ii) = ScalarBayesEstimators(M(1:l^n,1:n),wm(ii),...
                xmin,xmax,'method',method_opts,'estimator',estimator);
        end
        X = repmat(x{i}',[size(fBLS,1), 1, length(wm)]);
        Y = repmat(y{i}',[size(fBLS,1), 1, length(wm)]);
        fBLS = repmat(permute(fBLS,[1 3 2]),[1,size(X,2), 1]);
        WM = repmat(permute(wm(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
        WM_DRIFT = repmat(permute(wm_drift(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
        WY = repmat(permute(wy(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
        B = repmat(permute(b(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
        
        p_y_take_fBLS = (1./sqrt(2.*pi.*WY.^2.*fBLS.^2)) .* exp( -(Y - (fBLS+B)).^2./(2.*WY.^2.*fBLS.^2) );
        
        if n == 1
            p_m_take_x = (1./sqrt(2.*pi.*WM.^2.*X.^2)) .* ...
                exp( -squeeze(...
                (repmat(permute(M(1:l^n,1:n,:),[1 4 2 3]),[1 size(X,2) 1 1])...
                -repmat(permute(X,[1 2 4 3]),[1 1 n 1])).^2) ./...
                (2.*WM.^2.*X.^2) );
        elseif n == 2
            p_m1_take_x = (1./sqrt(2.*pi.*WM_DRIFT.^2.*X.^2)) .* ...
                exp( -squeeze(...
                (repmat(permute(M(1:l^n,1,:),[1 4 2 3]),[1 size(X,2) 1 1])...
                -repmat(permute(X,[1 2 4 3]),[1 1 n 1])).^2) ./...
                (2.*WM_DRIFT.^2.*X.^2) );
            p_m2_take_x = (1./sqrt(2.*pi.*WM.^2.*X.^2)) .*...
                exp( -squeeze(...
                (repmat(permute(M(1:l^n,2,:),[1 4 2 3]),[1 size(X,2) 1 1])...
                -repmat(permute(X,[1 2 4 3]),[1 1 n 1])).^2) ./...
                (2.*WM.^2.*X.^2) );
            p_m_take_x = p_m1_take_x.*p_m2_take_x;
        end
        
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
    estimator.type = 'SubOptMemBias';
    estimator.wm_drift = wm_drift;
    estimator.w_int = w_int;
    fBLS = nan(size(M(1:l^n,1:n),1),length(wm));
    for ii = 1:length(wm)
        fBLS(:,ii) = ScalarBayesEstimators(M(1:l^n,1:n),wm(ii),...
            xmin,xmax,'method',method_opts,'estimator',estimator);
    end
    X = repmat(x',[size(fBLS,1), 1, length(wm)]);
    Y = repmat(y',[size(fBLS,1), 1, length(wm)]);
    M = repmat(M,[1 1 length(wm)]);
    fBLS = repmat(permute(fBLS,[1 3 2]),[1,size(X,2), 1]);
    WM = repmat(permute(wm(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
    WM_DRIFT = repmat(permute(wm_drift(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
    WY = repmat(permute(wy(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
    B = repmat(permute(b(:),[2 3 1]),[size(fBLS,1), size(X,2), 1]);
    
    p_y_take_fBLS = (1./sqrt(2.*pi.*WY.^2.*fBLS.^2)) .* exp( -(Y - (fBLS+B)).^2./(2.*WY.^2.*fBLS.^2) );
    
    if n == 1
        p_m_take_x = (1./sqrt(2.*pi.*WM.^2.*X.^2)) .* ...
            exp( -squeeze(...
            (repmat(permute(M(1:l^n,1:n,:),[1 4 2 3]),[1 size(X,2) 1 1])...
            -repmat(permute(X,[1 2 4 3]),[1 1 n 1])).^2) ./...
            (2.*WM.^2.*X.^2) );
    elseif n == 2
        p_m1_take_x = (1./sqrt(2.*pi.*WM_DRIFT.^2.*X.^2)) .* ...
            exp( -squeeze(...
            (repmat(permute(M(1:l^n,1,:),[1 4 2 3]),[1 size(X,2) 1 1])...
            -repmat(permute(X,[1 2 4 3]),[1 1 n 1])).^2) ./...
            (2.*WM_DRIFT.^2.*X.^2) );
        p_m2_take_x = (1./sqrt(2.*pi.*WM.^2.*X.^2)) .*...
            exp( -squeeze(...
            (repmat(permute(M(1:l^n,2,:),[1 4 2 3]),[1 size(X,2) 1 1])...
            -repmat(permute(X,[1 2 4 3]),[1 1 n 1])).^2) ./...
            (2.*WM.^2.*X.^2) );
        p_m_take_x = p_m1_take_x.*p_m2_take_x;
    end
    
    integrand = p_y_take_fBLS.*p_m_take_x;
    
    logL = nan(length(wm),1);
    for ii = 1:length(wm)
        likelihood = (1-lambda(lapse(ii),Y,X)).*W(1:l^N)'*integrand(:,:,ii) + lambda(lapse(ii),Y,X).*uni(pmin,pmax,Y(:,:,ii));
        logL(ii,1) = -sum(log(likelihood),2);
    end
    
    
end



%% ModelEvidence functions
    function out = modelEvidenceFunction(ll,lfun)
        out = exp(-lfun + ll);
        out(isnan(out)) = 0;
        out = out(:);