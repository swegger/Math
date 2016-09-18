function [mu, V] = forwardBackwardPPS(x,model)
%% forwardBackwardPPS
%
%   [mu, V] = forwardBackwardPPS(y,model)
%
%   Computes an estimate of the posterior distribution for the GLMHS model
%   using a point-process adataptive filter (PPS). See Lawhern et. al. 2010
%   and Eden et. al. 2004).
%
%%

A = model.A;            % Transition matrix
Gamma = model.Gamma;    % Transition covariance

mu_c = model.mu_c;      % Mean firing rate
C = model.C;            % Linear dependence of spiking on hidden state
Sigma = model.Sigma;    % Neuron noise


hiddenDims = size(A,1);

%% Perform smoothing
% Initialize matrices
mu = zeros(hiddenDims,steps);
V = zeros(hiddenDims,hiddenDims,steps);

% Find estimate from kalman filter
[muhat, Vhat, K, mu_, V_] = kalmanFilter(x,A,C,Gamma,Sigma,mu0,V0);

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


%% Functions
function [mu, V] = forwardPass(x,mu0,V0,model)
    A = model.A;
    Gamma = model.Gamma;
    hiddenDims = size(A,1);
    
    % Initialize matrices
    mu_ = zeros(hiddenDims,steps);
    V_ = zeros(hiddenDims,hiddenDims,steps);
    mu = zeros(hiddenDims,steps);
    V = zeros(hiddenDims,hiddenDims,steps);
    
    % Iterate forward passes
    for k = 1:steps
        if k == 1
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
    
% Predictive distribution
function [mu, V] = predictive(mu0,V0,A,Gamma)
    mu = A*mu0;
    V = A*V0*A' + Gamma;