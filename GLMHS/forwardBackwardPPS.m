function [mu, V, Vj, mu_, V_, muhat, Vhat, J] = forwardBackwardPPS(x,model,varargin)
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
addParameter(Parser,'verbose',false)

parse(Parser,x,model,varargin{:})

x = Parser.Results.x;
model = Parser.Results.model;
steps = Parser.Results.steps;
mu0 = Parser.Results.mu0;
V0 = Parser.Results.V0;
verbose = Parser.Results.verbose;

if isnan(steps)
    steps = size(x,2);
end

if isnan(mu0)
    mu0 = randn(size(model.A,1),1);
end

if isnan(V0)
    V0 = eye(size(model.A,1));
end

%% Perform smoothing
% Generate forward estimates
if verbose
    disp('Performing forward passes...')
end
[muhat, Vhat, mu_, V_] = forwardPass(x,mu0,V0,model,steps);

% Perform backward recursion
if verbose
    disp('Performing backward passes...')
end
[mu, V, Vj, J] = backwardPass(muhat,mu_,Vhat,V_,model,steps);


%% Functions
% Foward Pass
function [mu, V, mu_, V_] = forwardPass(x,mu0,V0,model,steps)
    A = model.A;
    Gamma = model.Gamma;
    hiddenDims = size(A,1);
    C = model.C;
    mu_c = model.mu_c;
    
    
    % Initialize matrices
    mu_ = zeros(hiddenDims,steps);
    V_ = zeros(hiddenDims,hiddenDims,steps);
    mu = zeros(hiddenDims,steps);
    V = zeros(hiddenDims,hiddenDims,steps);
    
    % Iterate forward passes
    for k = 1:steps
        if k == 1
            % Compute predictive distribution;
            [mu_(:,k), V_(:,:,k)] = predictive(mu0,V0,A,Gamma);
            
            % Compute expected lambda
            l(:,k) = lambda(mu_(:,k),mu_c,C);
            for c = 1:size(x,1)
                temp(:,:,c) = C(c,:)'*l(c,k)*C(c,:);
                errs(:,c) = C(c,:)'*(x(c,k) - l(c,k));
            end
            Sigma = sum(temp,3);
            Err = sum(errs,2);

            % Update the predictive distribtuion using new measurements
            V(:,:,k) = inv( inv(V_(:,:,k)) + Sigma );
            mu(:,k) = mu_(:,k) + V(:,:,k)*Err;

        else
            % Compute predictive distribution;
            [mu_(:,k), V_(:,:,k)] = predictive(mu(:,k-1),V(:,:,k-1),A,Gamma);
            
            % Compute expected lambda
            l(:,k) = lambda(mu_(:,k),mu_c,C);
            for c = 1:size(x,1)
                temp(:,:,c) = C(c,:)'*l(c,k)*C(c,:);
                errs(:,c) = C(c,:)*(x(c,k) - l(c,k));
            end
            Sigma = sum(temp,3);
            Err = sum(errs,2);

            % Update the predictive distribtuion using new measurements
            V(:,:,k) = inv( inv(V_(:,:,k)) + Sigma );
            mu(:,k) = mu_(:,k) + V(:,:,k)*Err;

        end
        
        % Make sure covariance matrices are symetric
        V(:,:,k) = mean(cat(3,V(:,:,k),V(:,:,k)'),3);
        V_(:,:,k) = mean(cat(3,V_(:,:,k),V_(:,:,k)'),3);
        
        % Replace Infs w/ realmax
        V(isinf(V)) = realmax;
        V_(isinf(V_)) = realmax;
        
        % Make sure covariance matrices are positive definite
        [~,err] = cholcov(V(:,:,k));
        if err
            V(:,:,k) = nearestSPD(V(:,:,k));
        end
        [~,err] = cholcov(V_(:,:,k));
        if err
            V_(:,:,k) = nearestSPD(V_(:,:,k));
        end
        
    end

    
% Backward Pass
function [mu, V, Vj, J] = backwardPass(muhat,mu_,Vhat,V_,model,steps)

    A = model.A;
    mu = muhat;
    V = Vhat;
    
    for k = 1:steps-1
        % Update index
        indx = steps-k;

        % Find J
        J(:,:,indx) = Vhat(:,:,indx)*A'*inv(V_(:,:,indx+1));

        % Perform backwards pass
        mu(:,indx) = muhat(:,indx) + ...
            J(:,:,indx)*(mu(:,indx+1) - mu_(:,indx));
        V(:,:,indx) = Vhat(:,:,indx) + ...
            J(:,:,indx)*(V(:,:,indx+1) - V_(:,:,indx+1))*J(:,:,indx)';
    
        % Make sure covariance matrices are symetric
        V(:,:,k) = mean(cat(3,V(:,:,k),V(:,:,k)'),3);
        
        % Make sure covariance matrices are positive definite
        [~,err] = cholcov(V(:,:,k));
        if err
            V(:,:,k) = nearestSPD(V(:,:,k));
        end
    end
    
    Vj(:,:,steps) = V(:,:,steps)*J(:,:,end-1)';
    for k = 1:steps-2;
        % Update index
        indx = steps-k;
        
        % Find joint covariance
        Vj(:,:,indx) = V(:,:,indx)*J(:,:,indx-1)' + ...
            J(:,:,indx)*( Vj(:,:,indx+1) - A*Vhat(:,:,indx) )*J(:,:,indx-1)';
        
        % Make sure covariance matrices are symetric
%         Vj(:,:,k) = mean(cat(3,Vj(:,:,k),Vj(:,:,k)'),3);
        
        % Make sure covariance matrices are positive definite
%         [~,err] = cholcov(Vj(:,:,k));
%         if err
%             Vj(:,:,k) = nearestSPD(Vj(:,:,k));
%         end        
    end
    Vj(:,:,1) = V(:,:,1);
    
    
% Predictive distribution
function [mu, V] = predictive(mu0,V0,A,Gamma)
    mu = A*mu0;
    V = A*V0*A' + Gamma;