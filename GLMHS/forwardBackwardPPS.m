function [mu, V, Vj] = forwardBackwardPPS(x,model,varargin)
%% forwardBackwardPPS
%
%   [mu, V] = forwardBackwardPPS(y,model)
%
%   Computes an estimate of the posterior distribution for the GLMHS model
%   using a point-process adataptive filter (PPS). See Lawhern et. al. 2010
%   and Eden et. al. 2004).
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'x')
addRequired(Parser,'model')
addParameter(Parser,'steps',NaN)
addParameter(Parser,'mu0',NaN)
addParameter(Parser,'V0',NaN)

parse(Parser,x,model,varargin{:})

x = Parser.Results.x;
model = Parser.Results.model;
steps = Parser.Results.steps;
mu0 = Parser.Results.mu0;
V0 = Parser.Results.V0;

if isnan(steps)
    steps = size(x,2);
end

%% Perform smoothing
% Generate forward estimates
[muhat, Vhat, mu_, V_] = forwardPass(x,mu0,V0,model,steps);

% Perform backward recursion
[mu, V, Vj] = backwardPass(muhat,mu_,Vhat,V_,model,steps);


%% Functions
% Foward Pass
function [mu, V, mu_, V_] = forwardPass(x,mu0,V0,model,steps)
    A = model.A;
    Gamma = model.Gamma;
    hiddenDims = size(A,1);
    C = model.C;
    
    % Initialize matrices
    mu_ = zeros(hiddenDims,steps);
    V_ = zeros(hiddenDims,hiddenDims,steps);
    mu = zeros(hiddenDims,steps);
    V = zeros(hiddenDims,hiddenDims,steps);
    
    % Iterate forward passes
    for k = 1:steps
        if k == 1
            % Compute predictive distribution;
            [mu_(:,k), V_(:,:,k)] = predictive(mu0,V0,A,Gamma,1);
            
            % Compute expected lambda
            l(:,k) = lambda(mu_(:,k),mu_c,C);
            for c = 1:size(x,1)
                temp(:,:,c) = C(c,:)*l(c,k)*C(c,:)';
                errs(:,c) = C(c,:)*(x(c,k) - l(c,k));
            end
            Sigma = sum(temp,3);
            Err = sum(errs,2);

            % Update the predictive distribtuion using new measurements
            V(:,:,k) = inv( inv(V_(:,:,k)) + Sigma );
            mu(:,k) = mu_(:,k) + inv(V(:,:,k))*Err;

        else
            % Compute predictive distribution;
            [mu_(:,k), V_(:,:,k)] = predictive(mu(:,k-1),V(:,:,k-1),A,Gamma,1);
            
            % Compute expected lambda
            l(:,k) = lambda(mu_(:,k),mu_c,C);
            for c = 1:size(x,1)
                temp(:,:,c) = C(c,:)*l(c,k)*C(c,:)';
                errs(:,c) = C(c,:)*(x(c,k) - l(c,k));
            end
            Sigma = sum(temp,3);
            Err = sum(errs,2);

            % Update the predictive distribtuion using new measurements
            V(:,:,k) = inv( inv(V_(:,:,k)) + Sigma );
            mu(:,k) = mu_(:,k) + inv(V(:,:,k))*Err;

        end
    end

% Backward Pass
function [mu, V, Vj] = backwardPass(muhat,mu_,Vhat,V_,model,steps)

    A = model.A;
    
    for k = 1:steps-1
        % Update index
        indx = steps-k;

        % Find J
        J(:,:,indx) = Vhat(:,:,indx)*A'*V_(:,:,indx);

        % Perform backwards pass
        mu(:,indx) = muhat(:,indx) + ...
            J(:,:,indx)*(mu(:,indx+1) - mu_(:,indx));
        V(:,:,indx) = Vhat(:,:,indx) + ...
            J(:,:,indx)*(V(:,:,indx+1) - V_(:,:,indx))*J(:,:,indx)';

    end
    
    Vj(:,steps) = V(:,:,steps)*J(:,:,end-1)';
    for k = 1:steps-1;
        % Update index
        indx = steps-k;
        
        % Find joint covariance
        Vj(:,:,k) = V(:,:,k)*J(:,:,k-1)' + ...
            J(:,:,k)*( Vj(:,:,k+1) - A*Vhat(:,:,k) )*J(:,:,k-1)';
    end
    
    
% Predictive distribution
function [mu, V] = predictive(mu0,V0,A,Gamma)
    mu = A*mu0;
    V = A*V0*A' + Gamma;