function [mu, V, K, mu_, V_] = ...
    kalmanFilter(x,A,C,Gamma,Sigma,mu0,V0,varargin)
%% kalmanFilter
%
%   [mu, V] = kalmanFilter(x,A,Gamma,Sigma,mu0,V0)
%
%   Performs the Kalman filter algorithm to estimate the hidden state, z,
%   from observations, x. Model asssumes z evolves according to a linear
%   system of equations
%       z(k+1) = Az(k+1) + eta
%   where eta is distributed as a Gaussian with covariance Gamma.
%   Observations of the hidden state, x, are assumed to be a linear
%   function of the hidden state
%       x(k) = Cz(k) + eps
%   where eps is distributed as a Gaussian with covariance Sigma.
%
%   TODO:
%       (1) Assumes stationarity of measurement process and state
%       transition probabilities; future versions should allow for
%       nonstationarity.
%
%%

%% Defaults
a
%% Parse input
Parser = inputParser;

addRequired(Parser,'x')
addRequired(Parser,'A')
addRequired(Parser,'C')
addRequired(Parser,'Gamma')
addRequired(Parser,'Sigma')
addRequired(Parser,'mu0')
addRequired(Parser,'V0')
addParameter(Parser,'steps',Inf)

parse(Parser,x,A,C,Gamma,Sigma,mu0,V0,varargin{:})

x = Parser.Results.x;
A = Parser.Results.A;
C = Parser.Results.C;
Gamma = Parser.Results.Gamma;
Sigma = Parser.Results.Sigma;
mu0 = Parser.Results.mu0;
V0 = Parser.Results.V0;
steps = Parser.Results.steps;

if isinf(steps)
    steps = size(x,2);
end

hiddenDims = size(A,1);

%% Apply kalman filter for K steps from initial condition
mu_ = zeros(hiddenDims,steps);
V_ = zeros(hiddenDims,hiddenDims,steps);
mu = zeros(hiddenDims,steps);
V = zeros(hiddenDims,hiddenDims,steps);
for k = 1:steps
    if k == 1
        % Compute predictive distribution;
        [mu_(:,k), V_(:,:,k)] = kalmanPredictor(mu0,V0,A,Gamma,1);
        
        % Compute Kalman gain
        K(:,:,k) = kalmanGain(V_(:,:,k),C,Sigma);
        
        % Update the predictive distribtuion using new measurements
        mu(:,k) = mu_(:,k) + K(:,:,k)*(x(:,k) - C*mu_(:,k));
        V(:,:,k) = (eye(size(V0)) - K(:,:,k)*C)*V_(:,:,k);
        
    else
        % Compute predictive distribution;
        [mu_(:,k), V_(:,:,k)] = kalmanPredictor(mu(:,k-1),V(:,:,k-1),A,Gamma,1);
        
        % Compute Kalman gain
        K(:,:,k) = kalmanGain(V_(:,:,k),C,Sigma);
        
        % Update the predictive distribtuion using new measurements
        mu(:,k) = mu_(:,k) + K(:,:,k)*(x(:,k) - C*mu_(:,k));
        V(:,:,k) = (eye(size(V0)) - K(:,:,k)*C)*V_(:,:,k);
        
    end
end

%% Functions
% kalmanGain
function K = kalmanGain(V_,C,Sigma)
    K = V_*C'*inv(C*V_*C' + Sigma);