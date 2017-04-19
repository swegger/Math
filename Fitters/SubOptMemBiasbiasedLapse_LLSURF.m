function LLSURF = SubOptMemBiasbiasedLapse_LLSURF(Q,N,x,y,xmin,xmax,dx,pmin,pmax)
%%
%%

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

LLSURF = LLsurface(Q,N,x,y,xmin,xmax,dx,M,m,pmin,pmax);


%% FUNCTIONS


function logL = logLikelihoodQUAD(wm,wy,wm_drift,w_int,b,lapse,N,x,y,xmin,xmax,dx,M,m,pmin,pmax)
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
        p_y_take_fBLS( isnan(p_y_take_fBLS) ) = 0;
        
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
        
        integrand = p_y_take_fBLS.*p_m_take_x;
        
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
            -repmat(permute(X,[1 2 4 3]),[1 1 1 1])).^2) ./...
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

%% Likelihood surface function
        function LLSURF = LLsurface(Q,N,x,y,xmin,xmax,dx,M,m,pmin,pmax)
            
            % Cycle through parameter combinations
            parfor qi = 1:size(Q,1)
                disp(['Combination ' num2str(qi) ' of ' num2str (size(Q,1)])
                wm = Q(qi,1);
                wy = Q(qi,2);
                b = Q(qi,3);
                lapse = Q(qi,4);
                wm_drift = Q(qi,5);
                w_int = Q(qi,6);
                
                LLSURF(qi) = logLikelihoodQUAD(wm,wy,wm_drift,w_int,...
                    b,lapse,N,x,y,xmin,xmax,dx,M,m,pmin,pmax)
                
            end
            