%% TestNeuralNetworkEIF
%
%%

%% Variables

% General
dt = 0.0001;        % Size of time steps (s)
T = 1;              % Duration of experiment (s)
t = 0:dt:T-dt;         % Vector of time steps (s)
N = 2;            % Neurons
trials = 100;       % Number of trials

% Internal properties
urest = 0;          % Resting membrane potential
ureset = 0;         % Reset potential following spike
uT = 1;             % Threshold voltage for exponential nonlinearity
uS = 1;            % Sensitivity of exponential nonlinearity
R = 1;              % Input resistance
tau = 0.01;         % Membrane time-constant (C/gl: Capacitance per leak conductance)
thetaReset = 100;   % Membrane voltage at which a reset occurs
refractoryT = 0.005;% Refractory time of the neuron

% Cross-coupling
g = 0;                  % Gain term to control magnitude of cross-coupling
sigJ = g/sqrt(N);
J = sigJ*randn(N,N);    % Coupling matrix

% Input
I0 = 0.95;
sigI = 0.5;

%% Run simulation

% Initialize variables
u0 = rand(1,N);         % Inital membrane voltage
noise = sigI*randn(length(t),N,trials);
I = I0*ones(length(t),N,trials) + noise;
I2 = 0.99*ones(length(t),N,trials) + noise;

% Implement system
[du, u, spikes] = neuralNetworkEIF(u0,urest,ureset,uT,uS,R,tau,I,J,...
                  'T',T,'dt',dt,'thetaReset',thetaReset,...
                  'refractoryT',refractoryT,'trials',trials);
              
[du2, u2, spikes2] = neuralNetworkEIF(u0,urest,ureset,uT,uS,R,tau,I2,J,...
                     'T',T,'dt',dt,'thetaReset',thetaReset,...
                     'refractoryT',refractoryT,'trials',trials);
%% Analysis