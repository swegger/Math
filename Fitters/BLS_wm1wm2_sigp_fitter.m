function [wm, wy, wm_drift, b, lapse, sigp, llikelihood, lmodelEvidence] = ...
    BLS_wm1wm2_sigp_fitter(x,y,varargin)
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
ObsActflg = 0;
Boundsflg = 0;
InequalityBoundflg = 0;
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
        elseif strcmp(varargin{i},'CrossValidation')
            CrossValidationflg = 1;
            CrossValidationnum = i;
        elseif strcmp(varargin{i},'ModelEvidence')
            ModelEvidenceflg = 1;
            ModelEvidencenum = i;
        elseif strcmp(varargin{i},'ObsAct')
            ObsActflg = 1;
            ObsActnum = i;
        elseif strcmp(varargin{i},'Bounds')
            Boundsflg = 1;
            Boundsnum = i;
        elseif strcmp(varargin{i},'InequalityBound')
            InequalityBoundflg = 1;
            InequalityBoundnum = i;
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

if CrossValidationflg
    CrossValidation = varargin{CrossValidationnum+1};
else
    CrossValidation.Type = 'None';
end

if ModelEvidenceflg
    ModelEvidence = varargin{ModelEvidencenum+1};
    if ~isfield(ModelEvidence,'OpenMind')
        ModelEvidence.OpenMind = 0;
    end
else
    ModelEvidence.method = 'none';
end

if ObsActflg
    ObsAct = varargin{ObsActnum+1};
else
    ObsAct = 0;
end

if Boundsflg
    ub = varargin{Boundsnum+1}(1,:);
    lb = varargin{Boundsnum+1}(2,:);
else
    lb = [0.01 0 0.01 -Inf 0 0];
    ub = [1    1 1     Inf 1 Inf];
end

if InequalityBoundflg
    inequalityVar = varargin{InequalityBoundnum+1};
    if isstruct(inequalityVar)
        A = varargin{InequalityBoundnum+1}.A;
        C = varargin{InequalityBoundnum+1}.C;
    elseif islogical(inequalityVar) && inequalityVar
        A = [1 0 -1 0 0 0];       % if logical and true, assume wm - wm_drift < 0
        C = 0;
    end
else
    A = [];
    C = [];
end

if nargin < 3 
    error('Not enought input arguments!')
end

%% Set up cross validation
switch CrossValidation.Type
    case {'None','none'}
        % No cross validation, fit all the data
        xfit{1} = x;
        yfit{1} = y;
        xval{1} = x;
        yval{1} = y;
        
    case 'FitFraction'
        % Fit a fraction of the data and validate on the remaining
        if ~isfield(CrossValidation,'Iter')
            CrossValidation.Iter = 1;
        end
        for ii = 1:CrossValidation.Iter
            if iscell(x)
                for i = 1:length(x)
                    [xfit{ii}{i} indx] = datasample(x{i},ceil(CrossValidation.FitPercent*length(x{i})),'Replace',false);
                    yfit{ii}{i} = y{i}(indx);
                    valvec = zeros(length(x{i}),1);
                    valvec(indx) = 1;
                    valvec = ~valvec;
                    xval{ii} = x{i}(valvec);
                    yval{ii} = y{i}(valvec);
                end
                
            else
                [xfit{ii} indx] = datasample(x,ceil(CrossValidation.FitPercent*length(x)),'Replace',false);
                yfit{ii} = y(indx);
                valvec = zeros(length(x),1);
                valvec(indx) = 1;
                valvec = ~valvec;
                xval{ii} = x(valvec);
                yval{ii} = y(valvec);
            end
        end
        
    case 'LNOCV'
        % Leave N out crossvalidation: fit on all but N, validate on N,
        % repeat until all the data is processed
        if iscell(x)
            for i = 1:length(x)
                for ii = 1:ceil(length(x{i})/CrossValidation.N)
                    if (ii-1)*CrossValidation.N+CrossValidation.N < length(x{i})
                        indx = (ii-1)*CrossValidation.N+1:(ii-1)*CrossValidation.N+CrossValidation.N;
                        fitvec = true(1,length(x{i}));
                        fitvec(indx) = false; 
                        valvec = ~fitvec;
                        xfit{ii}{i} = x{i}( fitvec );
                        yfit{ii}{i} = y{i}( fitvec );
                        xval{ii}{i} = x{i}( valvec );
                        yval{ii}{i} = y{i}( valvec );
                    else
                        xfit{ii}{i} = x{i}( 1:(ii-1)*CrossValidation.N );
                        yfit{ii}{i} = y{i}( 1:(ii-1)*CrossValidation.N );
                        xval{ii}{i} = x{i}( (ii-1)*CrossValidation.N+1:end );
                        yval{ii}{i} = y{i}( (ii-1)*CrossValidation.N+1:end );
                    end
                    xfitsz{ii}(i,:) = size(xfit{ii}{i});
                end
            end
            iikeep = true(size(1:ceil(length(x{i})/CrossValidation.N)));
            for ii = 1:ceil(length(x{i})/CrossValidation.N)
                xfitu = unique(xfitsz{ii},'rows');
                if any(xfitu(:) == 0)
                    iikeep(ii) = false;
                end
            end
            xfit = xfit(iikeep);
            yfit = yfit(iikeep);
            xval = xval(iikeep);
            yval = yval(iikeep);
            
        else
            for ii = 1:ceil(length(x)/CrossValidation.N)
                if (ii-1)*CrossValidation.N+CrossValidation.N < length(x{i})
                    indx = (ii-1)*CrossValidation.N+1:(ii-1)*CrossValidation.N+CrossValidation.N;
                    fitvec = true(1,length(x));
                    fitvec(indx) = false;
                    valvec = ~fitvec;
                    xfit{ii} = x( fitvec );
                    yfit{ii} = y( fitvec );
                    xval{ii} = x( valvec );
                    yval{ii} = y( valvec );
                else
                    xfit{ii} = x( 1:(ii-1)*CrossValidation.N );
                    yfit{ii} = y( 1:(ii-1)*CrossValidation.N );
                    xval{ii} = x( (ii-1)*CrossValidation.N+1:end );
                    yval{ii} = y( (ii-1)*CrossValidation.N+1:end );
                end
            end
        end   
        
    otherwise
        error(['Cross validation type ' CrossValidation.Type ' not recognized!'])
end

% Fit according to minimizer
OPTIONS = optimset('Display','iter');

for ii = 1:length(xfit)
    disp(['Fit # ' num2str(ii) ' of ' num2str(length(xfit))])
    switch FitType
            
        case 'quad'
            % Use Simpson's quadrature
            m = 1:dx:2*xmax;
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
            
            minimizant = @(p)logLikelihoodQUAD(p(:,1),p(:,2),p(:,3),p(:,4),p(:,5),p(:,6),N,xfit{ii},yfit{ii},xmin,xmax,dx,M,m,LapseSupport(1),LapseSupport(2),ObsAct);
            validant = @(p)logLikelihoodQUAD(p(:,1),p(:,2),p(:,3),p(:,4),p(:,5),p(:,6),N,xval{ii},yval{ii},xmin,xmax,dx,M,m,LapseSupport(1),LapseSupport(2),ObsAct);
            
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
            
            minimizant = @(p)logLikelihoodQUADbatch(p(:,1),p(:,2),p(:,3),p(:,4),p(:,5),p(:,6),N,xfit{ii},yfit{ii},xmin,xmax,dx,M,m,LapseSupport(1),LapseSupport(2),batchsize,ObsAct);
            validant = @(p)logLikelihoodQUADbatch(p(:,1),p(:,2),p(:,3),p(:,4),p(:,5),p(:,6),N,xval{ii},yval{ii},xmin,xmax,dx,M,m,LapseSupport(1),LapseSupport(2),batchsize,ObsAct);
            
        case 'quad_nested'
            % Use Simpson's quadrature on batches of data, performing
            % integrations in nested loops
            m = 0:dx:2*xmax;
            l = length(m);
            
            minimizant = @(p)logLikelihoodQUADnested(p(:,1),p(:,2),p(:,3),p(:,4),p(:,5),p(:,6),N,xfit{ii},yfit{ii},xmin,xmax,dx,m,batchsize,ObsAct);
            validant = @(p)logLikelihoodQUADnested(p(:,1),p(:,2),p(:,3),p(:,4),p(:,5),p(:,6),N,xval{ii},yval{ii},xmin,xmax,dx,m,batchsize,ObsAct);
    end
    
    %minimizer = 'fminsearch(minimizant, [wM_ini wP_ini b_ini lapse_ini], OPTIONS);';
    minimizer = 'fmincon(minimizant, [wM_ini wP_ini wM_drift_ini b_ini lapse_ini sigp_ini], A, C, [], [], lb, ub, [], OPTIONS);';
    if ii == 1
        wM_ini = IC(1);
        wP_ini = IC(2);
        wM_drift_ini = IC(3);
        b_ini = IC(4);
        lapse_ini = IC(5);
        sigp_ini = IC(6);
    else
        wM_ini = wm(ii-1);
        wP_ini = wy(ii-1);
        b_ini = b(ii-1);
        lapse_ini = lapse(ii-1);
        wM_drift_ini = wm_drift(ii-1);
        sigp_ini = sigp(ii-1);
    end
   try
        [lparams, llike, exitflg, output, lambda, grad, hessian] = eval(minimizer);
   catch ME
       save('/om/user/swegger/BLS_wm1wm2_sigp_fitterERROR')
       rethrow(ME)
    end
    
    wm(ii) = lparams(1);
    wy(ii) = lparams(2);
    wm_drift(ii) = lparams(3);
    b(ii) = lparams(4);
    lapse(ii) = lparams(5);
    sigp(ii) = lparams(6);
    
    switch CrossValidation.Type
        case 'None'
            llikelihood(ii) = llike;
            
        otherwise
            llikelihood(ii) = validant(lparams);
            
    end
    
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

error('SubOptMemBias model not yet supported for FitType = "integral"')


function logL = logLikelihoodTRAPZ(wm,wy,b,N,m,x,y)
%% LOGLIKELIHOODTRAPZ
%
%   Calculates the loglikelihood of measurement scalar
%   varaibility (wm) and produciton scalar variability (wy), given the
%   observed sample (x) and production times (y)

error('BLS_wm1wm2_sigp model not yet supported for FitType = "trapz"')

function logL = logLikelihoodQUAD(wm,wy,wm_drift,b,lapse,sig,N,x,y,xmin,xmax,dx,M,m,pmin,pmax,ObsAct)
%% LOGLIKELIHOODQUAD
%
%   Calculates the log likelihood of scalar
%   varaibility (wm) and produciton scalar variability (wy), given the
%   observed sample (x) and production times (y) using Simpson's quadrature
%   method
%
%%

if size(N,2) > 2
    error('Support of models with larger than 2 measurements not yet supported when measurement noise is not equal')
end

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
        estimator.type = 'BLS_wm1wm2';
        estimator.wm_drift = wm_drift;
        estimator.ObsAct = ObsAct;
        estimator.wy = wy;
        
        f = nan(size(M(1:l^n,1:n),1),length(wm));
        for ii = 1:length(wm)
            f(:,ii) = ScalarBayesEstimators(M(1:l^n,1:n),wm(ii),...
                xmin,xmax,'method',method_opts,'estimator',estimator);
        end
        X = repmat(x{i}',[size(f,1), 1, length(wm)]);
        Y = repmat(y{i}',[size(f,1), 1, length(wm)]);
        f = repmat(permute(f,[1 3 2]),[1,size(X,2), 1]);
        WM = repmat(permute(wm(:),[2 3 1]),[size(f,1), size(X,2), 1]);
        WM_DRIFT = repmat(permute(wm_drift(:),[2 3 1]),[size(f,1), size(X,2), 1]);
        WY = repmat(permute(wy(:),[2 3 1]),[size(f,1), size(X,2), 1]);
        B = repmat(permute(b(:),[2 3 1]),[size(f,1), size(X,2), 1]);
        SIG = repmat(permute(sig(:),[2 3 1]),[size(f,1), size(X,2), 1]);
        
        p_y_take_f = (1./sqrt(2.*pi.*(WY.^2.*f.^2 +SIG.^2))) .* exp( -(Y - (f+B)).^2./(2.*(WY.^2.*f.^2 +SIG.^2)) );
        p_y_take_f( isnan(p_y_take_f) ) = 0;
        
        if n == 1
            p_m_take_x = (1./sqrt(2.*pi.*WM.^2.*X.^2)) .* ...
                exp( -squeeze(...
                (repmat(permute(M(1:l^n,1:n,:),[1 4 2 3]),[1 size(X,2) 1 1])...
                -repmat(permute(X,[1 2 4 3]),[1 1 n 1])).^2) ./...
                (2.*WM.^2.*X.^2) );
        elseif n == 2
%             p_m_take_x = (1./sqrt(2.*pi.*WM.^2.*X.^2)).^n .* exp( -squeeze(sum((repmat(permute(M(1:l^n,1:n,:),[1 4 2 3]),[1 size(X,2) 1 1])-repmat(permute(X,[1 2 4 3]),[1 1 n 1])).^2,3))./(2.*WM.^2.*X.^2) );
            p_m1_take_x = (1./sqrt(2.*pi.*WM_DRIFT.^2.*X.^2)) .* ...
                exp( -squeeze(...
                (repmat(permute(M(1:l^n,1,:),[1 4 2 3]),[1 size(X,2) 1 1])...
                -repmat(permute(X,[1 2 4 3]),[1 1 1 1])).^2) ./...
                (2.*WM_DRIFT.^2.*X.^2) );
            p_m2_take_x = (1./sqrt(2.*pi.*WM.^2.*X.^2)) .*...
                exp( -squeeze(...
                (repmat(permute(M(1:l^n,2,:),[1 4 2 3]),[1 size(X,2) 1 1])...
                -repmat(permute(X,[1 2 4 3]),[1 1 1 1])).^2) ./...
                (2.*WM.^2.*X.^2) );
            p_m_take_x = p_m1_take_x.*p_m2_take_x;
        end
        
        integrand = p_y_take_f.*p_m_take_x;
        
        for ii = 1:length(wm)
            likelihood = (1-lambda(lapse(ii),Y,X)).*W(1:l^n)'*integrand(:,:,ii) + lambda(lapse(ii),Y,X).*uni(pmin,pmax,Y(:,:,ii));
            logLi(i,ii) = -sum(log(likelihood),2);
        end
        
    end
    logL = permute(sum(logLi,1),[2 1]);

%     if any(isnan(logL(:)))
%         error('Returns NaN')
%     end
    
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
    estimator.type = 'BLS_wm1wm2';
    estimator.wm_drift = wm_drift;
    estimator.ObsAct = ObsAct;
    estimator.wy = wy;
        
    f = nan(size(M(1:l^N,1:N),1),length(wm));
    for ii = 1:length(wm)
        f(:,ii) = ScalarBayesEstimators(M(1:l^N,1:N),wm(ii),...
            xmin,xmax,'method',method_opts,'estimator',estimator);
    end
    X = repmat(x',[size(f,1), 1, length(wm)]);
    Y = repmat(y',[size(f,1), 1, length(wm)]);
    M = repmat(M,[1 1 length(wm)]);
    f = repmat(permute(f,[1 3 2]),[1,size(X,2), 1]);
    WM = repmat(permute(wm(:),[2 3 1]),[size(f,1), size(X,2), 1]);
    WM_DRIFT = repmat(permute(wm_drift(:),[2 3 1]),[size(f,1), size(X,2), 1]);
    WY = repmat(permute(wy(:),[2 3 1]),[size(f,1), size(X,2), 1]);
    B = repmat(permute(b(:),[2 3 1]),[size(f,1), size(X,2), 1]);
    SIG = repmat(permute(sig(:),[2 3 1]),[size(f,1), size(X,2), 1]);
    
    p_y_take_f = (1./sqrt(2.*pi.*(WY.^2.*f.^2 +SIG.^2))) .* exp( -(Y - (f+B)).^2./(2.*(WY.^2.*f.^2 +SIG.^2)) );
    
    if N == 1
        p_m_take_x = (1./sqrt(2.*pi.*WM.^2.*X.^2)) .* ...
            exp( -squeeze(...
            (repmat(permute(M(1:l^N,1:N,:),[1 4 2 3]),[1 size(X,2) 1 1])...
            -repmat(permute(X,[1 2 4 3]),[1 1 N 1])).^2) ./...
            (2.*WM.^2.*X.^2) );
    elseif N == 2
        p_m1_take_x = (1./sqrt(2.*pi.*WM_DRIFT.^2.*X.^2)) .* ...
            exp( -squeeze(...
            (repmat(permute(M(1:l^N,1,:),[1 4 2 3]),[1 size(X,2) 1 1])...
            -repmat(permute(X,[1 2 4 3]),[1 1 1 1])).^2) ./...
            (2.*WM_DRIFT.^2.*X.^2) );
        p_m2_take_x = (1./sqrt(2.*pi.*WM.^2.*X.^2)) .*...
            exp( -squeeze(...
            (repmat(permute(M(1:l^N,2,:),[1 4 2 3]),[1 size(X,2) 1 1])...
            -repmat(permute(X,[1 2 4 3]),[1 1 N 1])).^2) ./...
            (2.*WM.^2.*X.^2) );
        p_m_take_x = p_m1_take_x.*p_m2_take_x;
    end
    
    integrand = p_y_take_f.*p_m_take_x;
    
    logL = nan(length(wm),1);
    for ii = 1:length(wm)
        likelihood = (1-lambda(lapse(ii),Y,X)).*W(1:l^N)'*integrand(:,:,ii) + lambda(lapse(ii),Y,X).*uni(pmin,pmax,Y(:,:,ii));
        logL(ii,1) = -sum(log(likelihood),2);
    end
    
    
end


function logL = logLikelihoodQUADbatch(wm,wy,b,lapse,N,x,y,xmin,xmax,dx,M,m,pmin,pmax,batchsize,ObsAct)
%% LOGLIKELIHOODQUADbatch
%
%   Calculates the log likelihood of scalar
%   varaibility (wm) and produciton scalar variability (wy), given the
%   observed sample (x) and production times (y) using Simpson's quadrature
%   method
%
%%

error('Lapse model not yet supported for FitType = "quad_batch"')

function logL = logLikelihoodQUADnested(wm,wy,b,N,x,y,xmin,xmax,dx,m,batchsize,ObsAct)


error('Lapse model not yet supported for FitType = "quad_nested"')

%% ModelEvidence functions
    function out = modelEvidenceFunction(ll,lfun)
        out = exp(-lfun + ll);
        out(isnan(out)) = 0;
        out = out(:);