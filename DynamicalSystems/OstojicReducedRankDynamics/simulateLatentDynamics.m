function [dkappas, kappas, sigmaTildes, deltas, vs] = simulateLatentDynamics(varargin)
%%
%
%
%
%%

%% Defaults
plotOpts_default.On = false;

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'tau',20)
addParameter(Parser,'t',0:2000)
addParameter(Parser,'us',NaN)
addParameter(Parser,'kappas0',NaN)
addParameter(Parser,'overlaps',NaN)
addParameter(Parser,'sigmas',NaN)
addParameter(Parser,'forceLinear',false)
addParameter(Parser,'plotOpts',plotOpts_default)

parse(Parser,varargin{:})

tau = Parser.Results.tau;
t = Parser.Results.t;
us = Parser.Results.us;
kappas0 = Parser.Results.kappas0;
overlaps = Parser.Results.overlaps;
sigmas = Parser.Results.sigmas;
forceLinear = Parser.Results.forceLinear;
plotOpts = Parser.Results.plotOpts;

%% Check input, generate if necessary
if any(isnan(us(:)))
    us = zeros(2,size(t,2));
    us(1,t>=500) = 0.1;
    
    us(2,t>=1000) = -1;
end
M = size(us,1);

%% Initialize kappas
if any(isnan(kappas0))
    kappas0(1,:) = 0;
    kappas0(2,:) = 0;
end
N = size(kappas0,1);

%% Initialize overlaps and sigmas;
if any(isnan(overlaps(:)))
    overlaps = zeros(N,N+M);
    overlaps(1,1) = 2.7;
    overlaps(2,1) = 2.0;
    overlaps(2,2) = 0.001;
    overlaps(1,N+1) = 0.2;
    overlaps(2,N+2) = 0;
    overlaps(2,N+1) = -1.2;
end

if any(isnan(sigmas))
    sigmas = 1.5*ones(N+M,1); %randn(N+M,1);
end

%% Run simulation
[dkappas, kappas, sigmaTildes, deltas, vs] = latenDynamicsRR(tau,us,overlaps,kappas0,sigmas,'forceLinear',forceLinear);

%% Plot results
if plotOpts.On
    figure;
    subplot(4,1,1)
    plot(t,vs)
    xlabel('Time (ms)')
    ylabel('Inputs')
    
    subplot(4,1,2)
    plot(t,dkappas')
    xlabel('Time (ms)')
    ylabel('d\kappa/dt')
    
    subplot(4,1,3)
    plot(t,kappas')
    xlabel('Time (ms)')
    ylabel('\kappa')
    
    subplot(4,1,4)
    plot(t,deltas)
    xlabel('Time (ms)')
    ylabel('\Delta')
end
