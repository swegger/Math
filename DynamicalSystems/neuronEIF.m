function [du, u, spikes] = neuronEIF(u0,urest,ureset,uT,uS,R,tau,I,varargin)
%% neuronEIF
%
%   [du, u, spikes] = neuronEIF(u0,urest,ureset,uT,uS,R,tau,I)
%   
%   Simulates the exponential-integrate-and-fire neuron model according to
%   the parameters, usering the Euler method:
%       u0 - initial conditions of the neuron voltage
%       urest - resting voltage of the neuron (leak potential)
%       ureset - reset voltage following a spike
%       uT - Rheobase threshold
%       uS - "Sharpness" of spike initiation process
%       R - input resistence
%       tau - membrane time constant
%       I - input to the neuron. Columns are different trials, rows are
%           input over time. If size(I,1) == 1, assume constant input.
%
%   ... = neuronEIF(...,'T',T)
%       Allows user to specify the total duration of the simulation. If
%       size(I,1) == 1, simulation assumes constant input for T/dt time
%       steps. If size(I,1) ~= 1 & size(I,1) ~= T/dt, error will occur.
%
%   ... = neuronEIF(...,'dt',dt)
%       Allows user to specify the duration of a timestep. Defaults to 1.
%
%   ... = neuronEIF(...,'thetaReset',thetaReset)
%       Allows user to specify the voltage at which a spike occurs. Should
%       be >> uT + uS.
%
%   ... = neuronEIF(...,'refractoryT',refractoryT)
%       Allows user to specify the abolute refractory time of the neuron.
%
%   Refs:
%       Fourcaud-Trocme N, et. al., 2003
%       Gerstner W, Kistler WM, Naud R, Paninski L. 2014
% v1 swe
%%

%% Defaults

%% Parse inputs
Parser = inputParser;
addRequired(Parser,'u0');
addRequired(Parser,'urest');
addRequired(Parser,'ureset');
addRequired(Parser,'uT');
addRequired(Parser,'uS');
addRequired(Parser,'R');
addRequired(Parser,'tau');
addRequired(Parser,'I');
addParameter(Parser,'T',[])
addParameter(Parser,'dt',1)
addParameter(Parser,'thetaReset',[])
addParameter(Parser,'refractoryT',[])

parse(Parser,u0,urest,ureset,uT,uS,R,tau,I,varargin{:})

u0 = Parser.Results.u0;
urest = Parser.Results.urest;
ureset = Parser.Results.ureset;
uT = Parser.Results.uT;
uS = Parser.Results.uS;
R = Parser.Results.R;
tau = Parser.Results.tau;
I = Parser.Results.I;
T = Parser.Results.T;
dt = Parser.Results.dt;
thetaReset = Parser.Results.thetaReset;
refractoryT = Parser.Results.refractoryT;

% Set initial values
if numel(u0) == 1
    u0 = repmat(u0,[1 size(I,2)]);
elseif numel(u0) ~= size(I,2)
    error('Length of initial conditions (u0) should match number of trials (size(I,2)!')
end

% Set the input/duration
if isempty(T)
    T = size(I,1)/dt;
elseif size(I,1) == 1
    I = repmat(I,[T/dt 1]);
elseif size(I,1) ~= T/dt
    error('Size of input must match duration of simulation divided by time step magnitude (dt)')
end

% Set thetaReset, if not provided
if isempty(thetaReset)
    thetaReset = 10 * uT+uS;
end

% Set the absolute refractory period of the neuron, if not provided
if isempty(refractoryT)
    refractoryT = 5*dt;
end

%% Initialize variables
du = nan(size(I));
u = nan(size(I));
spikes = nan(size(I));

u(1,:) = u0;

%% Simulate the model

for i = 2:T/dt
    % Determine the change in voltage
    du(i,:) = dt/tau * (...
        - (u(i-1,:) - urest)...                 % Leak current
        + uS*exp( (u(i-1,:) - uT)/uS ) ...      % Exponential nonlinearity
        + R*I(i,:) ...                          % Input current
        );
    
    % Update voltage
    u(i,:) = u(i-1,:) + du(i,:);    
        
    % Detect spikes and reset voltage
    spikes(i,:) = u(i,:) >= thetaReset;
    if any(spikes(i,:))
        disp(u(i-1:i,:))
    end
    
    % Instantiate voltage reset and absolute refractory period
    if i > floor(refractoryT/dt)
        u(i,(sum(spikes(i-floor(refractoryT/dt):i,:),1) > 0)) = ureset;
    else
        u(i,(sum(spikes(1:i,:),1) > 0)) = ureset;
    end
end