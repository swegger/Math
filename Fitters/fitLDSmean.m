function [x0, A, B, xhat, yhat] = fitLDSmean(y,N,varargin)
%% FitLDSmean
%
%   [x0, A, B, xhat, yhat] = FitLDSmean(y,N);
%      Fits a linear dynamical system to the data of the form
%           x(i+1) = x(i) + A*x(i);
%           y(i) = B*x(i);
%
%      Where x(i) is an N-dimensional vector of hidden states at step i and
%      y(i) is an M-dimensional vector of obervables. A is NxN and B is
%      MxN.
%
%      Theta is mapped as follow
%           x(1,:) = Theta(1:N)
%           A = [Theta(N+1) Theta(N+2) ... Theta(N+1+N)
%                Theta(N+1+N+1)...;
%                ...
%                Theta(N+1+N^2-N+1) ... Theta(N+1+N^2)]'
%           B = [Theta(N+1+N^2+1) Theta(N+1+N^2+2) ... Theta(N+1+N^2+M)
%                Theta(N+1+N^2+1+M+1) ...
%                ...
%                Theta(N+1+N^2+N*M-M+1) ... Theta(N+1+N^2+N*M)]'
%           
%
% swe
%%

%% Defaults
optionsDefault = optimset('Display','iter');

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'y')
addRequired(Parser,'N')
addParameter(Parser,'Theta0',NaN)
addParameter(Parser,'lowerBound',NaN)
addParameter(Parser,'upperBound',NaN)
addParameter(Parser,'constraintA',[])
addParameter(Parser,'constraintb',[])
addParameter(Parser,'options',optionsDefault)

parse(Parser,y,N,varargin{:});

y = Parser.Results.y;
N = Parser.Results.N;
Theta0 = Parser.Results.Theta0;
lowerBound = Parser.Results.lowerBound;
upperBound = Parser.Results.upperBound;
constraintA = Parser.Results.constraintA;
constraintb = Parser.Results.constraintb;
options = Parser.Results.options;
M = size(y,2);
steps = size(y,1);

if isnan(Theta0)
    Theta0 = -rand(N + N^2+N*M,1);      % N initial conditions, NxN interactions, and NxM observation weights
end

if isnan(lowerBound)
    lowerBound = -inf(size(Theta0));
end
if isnan(upperBound)
    upperBound = inf(size(Theta0));
end

%% Set up obective function
minimizant = @(p)( minFun(p,y,N,M,steps) );

%% Find parameters that minimize squared errors
% Minimize errors
Theta = fmincon(minimizant,Theta0,constraintA,constraintb,[],[],lowerBound,upperBound,[],options);

%% Find the predicted values of the hidden and observable states
[yhat, xhat] = LDSsim(Theta,N,M,steps);

%% Unpack parameter vector
x0 = Theta(1:N);        % Initial conditions
A = reshape(Theta(N+1:N+1+N^2-1),[N N]);      % Transition matrix
B = reshape(Theta(N+1+N^2:end),[M N]);    % Observation matrix

%% Functions
% For simulating the system
function [y, x] = LDSsim(Theta,N,M,steps)
% Set up parameters
x(1,:) = Theta(1:N);        % Initial conditions
A = reshape(Theta(N+1:N+1+N^2-1),[N N]);      % Transition matrix
B = reshape(Theta(N+1+N^2:end),[M N]);    % Observation matrix

y(1,:) = B*x(1,:)';

% Run the simulation
for i = 2:steps
    x(i,:) = x(i-1,:)' + A*x(i-1,:)';
    y(i,:) = B*x(i,:)';
end

% Minimization function
function err = minFun(Theta,y,N,M,steps)
yhat = LDSsim(Theta,N,M,steps);
err = sum( (y(:)-yhat(:)).^2 );