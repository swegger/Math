function [model, LL, z, V, Vj, Vhat, V_, muhat, mu_] = ...
    identifyGLMHS_EM(x,model0,varargin)
%% identifyGLMHS_EM
%
%   [model, LL, z] = idenfityGLMHS_EM(x,model0)
%
%   Finds the maximum likelihood model parameters for the GLMHS model using
%   expectation-maximization (EM) algorithm. See Lawhern et. al. 2010.
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'x')             % Observed spiking activity
addRequired(Parser,'model0')        % Initial model parameters
addParameter(Parser,'tolerance',1)    % Convergence tolerance
addParameter(Parser,'samples',10)   % Random samples to use per time step in Monte Carlo integration of likelihood function
addParameter(Parser,'initialParams',NaN)    % Initial parameter values to 
addParameter(Parser,'mu0',NaN);             % Initial condition of hidden variables
addParameter(Parser,'V0',NaN);              % A priori covariance in hidden variables
addParameter(Parser,'verbose',false)        % Display progress

parse(Parser,x,model0,varargin{:})

x = Parser.Results.x;
model0 = Parser.Results.model0;
tolerance = Parser.Results.tolerance;
samples = Parser.Results.samples;
mu0 = Parser.Results.mu0;
V0 = Parser.Results.V0;
verbose = Parser.Results.verbose;

if isnan(mu0)
    mu0 = randn(size(model0.A,1),size(x,2));
end

if isnan(V0)
    V0 = eye(size(model0.A,1));
    V0 = repmat(V0,[1,1,size(x,2)]);
end

%% EM
% Initialize values
model(1) = model0;
indx = 0;
LL(1) = -realmax;
difference = Inf;

% Perform EM
while difference > tolerance
    indx = indx+1;
    
    % Perform expectation step using forward-backward algorithm and
    % Gaussian approximation (point-process filter; see Edan et. al. 2005)
    if verbose
        disp('Peforming expectation step...')
    end
    for triali = 1:size(x,3)
        [z(:,:,triali), V(:,:,:,triali), Vj(:,:,:,triali), mu_(:,:,triali),...
            V_(:,:,:,triali), muhat(:,:,triali), Vhat(:,:,:,triali),...
            J(:,:,:,triali)] = forwardBackwardPPS(x(:,:,triali),model(indx),...
            'mu0',mu0(:,:,triali),'V0',V0(:,:,triali));
    end
    
    % Perform maximization step
    if verbose
        disp('Peforming maximization step...')
    end
    model(indx+1).A = ID_A(z,V,Vj);                             % Identify linear transition matrix
    model(indx+1).Gamma = ID_Gamma(z,model(indx+1).A,V_,Vhat,Vj,V,J);       % Identify transition noise %% TODO: CHECK IF RIGHT %%
%    model(indx+1).Gamma = model(indx).Gamma;   

    initialParams(1,:) = model(indx).mu_c;
    initialParams(2:size(model(indx).C,2)+1,:) = model(indx).C';
    [model(indx+1).mu_c, model(indx+1).C] = maxE2(x,z,V,initialParams,verbose);    % Maximization step for neural parameters
    
    % Calculate log likelihood of parameters, given data
    LL(indx+1) = logLikelihood(x,mu_,V_,samples,model(indx+1).mu_c,model(indx+1).C);
    
    % Check for convergence
    difference = abs(LL(indx+1)-LL(indx));
    
    % Display output
    disp([num2str(indx) '    ' num2str(LL(indx+1)) '    ' num2str(difference)])
end

%% Functions

% ID_A
function A = ID_A(z,V,Vj)
    indx = 0;
    for triali = 1:size(z,3)
        for k = 1:size(z,2)-1
            indx = indx+1;
            numerator(:,:,indx) = Vj(:,:,k+1,triali) + z(:,k+1,triali)*z(:,k,triali)';
            denominator(:,:,indx) = V(:,:,k,triali) + z(:,k,triali)*z(:,k,triali)';
        end
    end
    A = sum(numerator,3) * inv(sum(denominator,3));


% ID_Gamma %% TODO: CHECK IF RIGHT %%
function Gamma = ID_Gamma(z,A,V_,Vhat,Vj,V,J)
    indx = 0;
    for triali = 1:size(z,3)
        for k = 2:size(z,2)
           indx = indx+1;
           expZn1Zn2(:,:,indx) = J(:,:,k-1,triali)*V(:,:,k,triali) + z(:,k,triali)*z(:,k-1,triali)';
    %         Gammas(:,:,k) = V_(:,:,k) - A*Vhat(:,:,k-1)*A';       % Estimates transition noise from changes in estimate covariance
    %         expZn1Zn2(:,:,k) = Vj(:,:,k) + z(:,k)*z(:,k-1)';       % Bishop, 2006 pg. 643
            expZn2Zn1(:,:,indx) = Vj(:,:,k,triali)' + z(:,k-1,triali)*z(:,k,triali)';      % Bishop, 2006 pg. 643
            expZn1Zn1(:,:,indx) = V(:,:,k-1,triali) + z(:,k-1,triali)*z(:,k-1,triali)';    % Bishop, 2006 pg. 643
            expZn2Zn2(:,:,indx) = V(:,:,k,triali) + z(:,k,triali)*z(:,k,triali)';          % Bishop, 2006 pg. 643
            Gammas(:,:,indx) = expZn2Zn2(:,:,indx) - A*expZn2Zn1(:,:,indx)...
                 - expZn1Zn2(:,:,indx)*A' + A*expZn1Zn1(:,:,indx)*A'; % Bishop, 2006 pg. 643
    %        Gammas(:,:,k) = expZn2Zn2 - A*expZn1Zn2' - expZn1Zn2*A + A*expZn1Zn1*A'; % Bishop, 2006 pg. 643
        end
    end
%     Gamma = mean(Gammas,3);     % Assumes stationarity of transition noise
    Gamma = 1/((size(z,2)-1)*size(z,3)) * sum(Gammas,3);    % Bishop, 2006 pg. 643
    Gamma = (Gamma + Gamma')/2;
    [~,err] = cholcov(Gamma);
    if err
        Gamma = nearestSPD(Gamma);
    end

    
% maxE2
function [mu_c, C] = maxE2(x,z,V,initialParams,verbose)
    
    minimizant = @(P)(E2(x,z,V,initialParams,P));
    initCond = initialParams';
    %initCond = initCond(:) + randn(size(initCond(:)))/100;
    if verbose
        OPTIONS = optimset('Display','iter','TolFun',20,'TolX',1000,'PlotFcns',@optimplotfval);
    else
        OPTIONS = optimset('Display','off','TolFun',5,'TolX',1000);
    end
    P = fminsearch(minimizant,initCond,OPTIONS);
    mu_c = P(1:size(x,1));
    mu_c = mu_c(:);
    C = reshape(P(size(x,1)+1:end),[size(x,1),size(z,1)]);

%     tol = 1;
%     % Perform Newton-Rahpson to find roots
%     for c = 1:size(x,1)     % unit loop
%         
%         % Initialize values
%         indx = 0;
%         Params(:,1) = initialParams(:,c);
%         d = Inf;
%         
%         while d > tol
%             indx = indx+1;
%             
%             mu_cTemp = Params(1,indx);
%             l = Params(2:end,indx);
%             
%             for k = 2:size(x,2) % time step loop
%                 S1 = mu_cTemp;
%                 S2 = l'*z(:,k) + 1/2*l'*V(:,:,k)*l;
%                 S3 = z(:,k) + V(:,:,k)*l;
%                 temp = zeros(size(V(:,:,k),1)+1,size(V(:,:,k),2)+1);
%                 temp(end-size(V(:,:,k),1)+1:end,end-size(V(:,:,k),2)+1:end) = V(:,:,k);
%                 Js(:,:,k) = exp(S1+S2)*( [1; S3]*[1; S3]' + temp );
%                 fs(:,k) = -exp(S1)*exp(S2)*[1; S3] - x(c,k)*[1; z(:,k)];
%             end
%             J = -sum(Js,3); 
%             f = sum(fs,2);
%             Params(:,indx+1) = Params(:,indx) - J\f;
%             
%             disp(f)
%             d = sum( (Params(:,indx+1) - Params(:,indx)).^2 );
%         end
%         mu_c(c,1) = Params(1,indx+1);
%         C(c,:) = Params(2:end,indx+1);
%         
%     end
    
    
% logLikelihood
function LL = logLikelihood(x,mu_,V_,samples,mu_c,C)

    % For each time point, estimate the log likelihood of x, given the
    % model parameters
    for triali = 1:size(x,3)
        for k = 1:size(x,2)
            Generate samples of hidden varialbe
            q_samp = mvnrnd(mu_(:,k,triali),V_(:,:,k,triali),samples);
            for sampi = 1:samples
                for c = 1:size(x,1)
                    likes(sampi,c) = 1/samples * (lambda(q_samp(sampi,:)',mu_c(c),C(c,:)))^x(c,k,triali)*...
                        exp(-lambda(q_samp(sampi,:)',mu_c(c),C(c,:))) / factorial(x(c,k,triali)) ;
                end
            end
            logLikes(k,triali) = log( 1/samples * sum ( prod(likes,2) , 1) );

        end
    end
    
    
    % Combine likelihoods across time and trials
    LL = sum( logLikes(:) );
    

% E2
function out = E2(x,z,V,P0,P)
    mu_c0 = P0(1,:)';
    C0 = P0(2:end,:)';
    ps = reshape(P,[size(x,1),size(z,1)+1]);
    mu_c = ps(:,1);
    C = ps(:,2:end);

    for triali = 1:size(x,3)
        for k = 1:size(x,2)
            for c = 1:size(x,1)
                for i = 1:10
                    temp = nearestSPD(V(:,:,k,triali));
                    zi = mvnrnd(z(:,k,triali),V(:,:,k,triali));
                    l = lambda(zi(:),mu_c(c),C(c,:));
                    monteCarloE(i) = -l + x(c,k,triali)*log(l);
                end
                e2(c,k,triali) = mean(monteCarloE);
            end
        end
    end

    out = -sum(e2(:));