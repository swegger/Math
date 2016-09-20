function [model, LL] = identifyGLMHS_EM(x,model0,varargin)
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

parse(Parser,x,model0,varargin{:})

x = Parser.Results.x;
model0 = Parser.Results.model0;
tolerance = Parser.Results.tolerance;
samples = Parser.Results.samples;

%% EM
% Initialize values
model(1) = model0;
indx = 0;
LL(1) = 0;
difference = Inf;

% Perform EM
while difference > tolerance
    indx = indx+1;
    
    % Perform expectation step using forward-backward algorithm and
    % Gaussian approximation (point-process filter; see Edan et. al. 2005)
    [z, V, Vj, mu_, V_, muhat, Vhat] = forwardBackwardPPS(x,model(indx));
    
    % Perform maximization step
    model(indx+1).A = ID_A(z,V,Vj);                             % Identify linear transition matrix
    model(indx+1).Gamma = ID_Gamma(z,model(indx+1).A,V_,Vhat);       % Identify transition noise
    
    initialParams(1,:) = model(indx).mu_c;
    initialParams(2:size(model(indx).C,2)+1,:) = model(indx).C';
    [model(indx+1).mu_c, model(indx+1).C] = maxE2(x,z,V,initialParams);    % Maximization step for neural parameters
    
    % Calculate log likelihood of parameters, given data
    LL(indx+1) = logLikelihood(x,mu_,V_,samples,model(indx+1).mu_c,model(indx+1).C);
    
    % Check for convergence
    difference = LL(indx+1)-LL(indx);
    
    % Display output
    disp([num2str(indx) '    ' num2str(LL(indx+1)) '    ' num2str(difference)])
end

%% Functions

% ID_A
function A = ID_A(z,V,Vj)
    for k = 1:size(z,2)-1
        numerator(:,:,k) = Vj(:,:,k+1) + z(:,k+1)*z(:,k)';
        denominator(:,:,k) = V(:,:,k) + z(:,k)*z(:,k)';
    end
    A = sum(numerator,3) * inv(sum(denominator,3));


% ID_Gamma %% TODO: CHECK IF RIGHT %%
function Gamma = ID_Gamma(z,A,V_,Vhat)
    for k = 2:size(z,2)
        Gammas(:,:,k) = V_(:,:,k) - A*Vhat(:,:,k-1)*A';       % Estimates transition noise from changes in estimate covariance
    end
    Gamma = mean(Gammas,3);     % Assumes stationarity of transition noise

    
% maxE2
function [mu_c, C] = maxE2(x,z,V,intialParams)

    tol = 1;
    % Perform Newton-Rahpson to find roots
    for c = 1:size(x,1)     % unit loop
        
        % Initialize values
        indx = 0;
        Params(:,1) = intialParams(:,c);
        d = Inf;
        
        while d > tol
            indx = indx+1;
            
            mu_c = Params(1,indx);
            l = Params(2:end,indx);
            
            for k = 2:size(x,2) % time step loop
                S1 = mu_c;
                S2 = l'*z(:,k) + 1/2*l'*V(:,:,k)*l;
                S3 = z(:,k) + V(:,:,k)*l;
                temp = zeros(size(V(:,:,k),1)+1,size(V(:,:,k),2)+1);
                temp(end-size(V(:,:,k),1)+1:end,end-size(V(:,:,k),2)+1:end) = V(:,:,k);
                Js(:,:,k) = exp(S1+S2)*( [1; S3]*[1; S3]' + temp );
                fs(:,k) = -exp(S1)*exp(S2)*[1; S3] - x(c,k)*[1; z(:,k)];
            end
            J = -sum(Js,3);
            f = sum(fs,2);
            Params(:,indx+1) = Params(:,indx) - inv(J)*f;
            
            d = sum( (Params(:,indx+1) - Params(:,indx)).^2 );
        end
        mu_c(c,1) = Params(1,indx+1);
        C(c,:) = Params(2:end,indx+1);
        
    end
    
    
% logLikelihood
function LL = logLikelihood(x,mu_,V_,samples,mu_c,C)

    % For each time point, estimate the log likelihood of x, given the
    % model parameters
    for k = 1:size(x,2)
        % Generate samples of hidden varialbe
        q_samp = mvnrnd(mu_(:,k),V_(:,:,k),samples);
        for sampi = 1:samples
            for c = 1:size(x,1)
                likes(sampi,c) = 1/samples * (lambda(q_samp(sampi,:)',mu_c(c),C(c,:)))^x(c,k)*...
                    exp(-lambda(q_samp(sampi,:)',mu_c(c),C(c,:))) / factorial(x(c,k)) ;
            end
        end
        logLikes(k) = log( 1/samples * sum ( prod(likes,2) , 1) );
        
    end
    
    
    % Combine likelihoods across time
    LL = sum( logLikes );