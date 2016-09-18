function [mu, V, muhat, Vhat, mu_, V_, J, K] = ...
    kalmanSmoother(x,A,C,Gamma,Sigma,mu0,V0,varargin)
%% kalmanSmoother
%
%   [mu, V] = kalmanSmoother(x,A,C,Gamma,Sigma,mu0,V0)
%
%   Computes the optimal estimate of the hidden state at each time point,
%   z(k), using the full set of observations x.
%
%%

%% Defaults

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
addParameter(Parser,'algorithm','ForwardBackward')

parse(Parser,x,A,C,Gamma,Sigma,mu0,V0,varargin{:})

x = Parser.Results.x;
A = Parser.Results.A;
C = Parser.Results.C;
Gamma = Parser.Results.Gamma;
Sigma = Parser.Results.Sigma;
mu0 = Parser.Results.mu0;
V0 = Parser.Results.V0;
steps = Parser.Results.steps;
algorithm = Parser.Results.algorithm;

if isinf(steps)
    steps = size(x,2);
end

hiddenDims = size(A,1);

%% Perform smoothing
switch algorithm
    case {'ForwardBackward','Rauch-Tung-Striebel','RTS'}
        % Initialize matrices
        mu = zeros(hiddenDims,steps);
        V = zeros(hiddenDims,hiddenDims,steps);
        J = zeros(hiddenDims,hiddenDims,steps);
        
        % Find estimate from kalman filter
        [muhat, Vhat, K, mu_, V_] = kalmanFilter(x,A,C,Gamma,Sigma,mu0,V0);
        
        mu(:,steps) = muhat(:,steps);
        V(:,steps) = Vhat(:,steps);
        for k = 1:steps-1
            % Update index
            indx = steps-k;
            
            % Find J
            J(:,:,indx) = Vhat(:,:,indx)*A'*inv(V_(:,:,indx));
            
            % Perform backwards pass
            mu(:,indx) = muhat(:,indx) + ...
                J(:,:,indx)*(mu(:,indx+1) - mu_(:,indx));
            V(:,:,indx) = Vhat(:,:,indx) + ...
                J(:,:,indx)*(V(:,:,indx+1) - V_(:,:,indx))*J(:,:,indx)';
            
        end

    case 'BlockTridiagonal' % See Paninskit et. al. 2010
        error('Not yet implemented')
        D = inv(Gamma) + A'*inv(Gamma)*A + C'*inv(Sigma)*C;
        R = -inv(Gamma)*A;
        
        H = nan;
        
    otherwise
        error('Kalman smoothing algorithm not recognized!')
        
end
