%% TestFitDialationFactor
%
%   Script generating spikes or analog data from a process that scales with
%   time and a process that does not scale with time.
%
%%

%% Variables
% General
durs = 600:200:1400;                        % Duration of trials
trials = 100;                               % Trials per neuron
t = 0:max(durs);                                    % Time vector

% Dynamics
tau = 200;                                  % Time constant of unscaled data.
taus = tau*durs(1)/max(durs) + (tau - tau*durs(1)/max(durs))./(1 + 10.^(-0.005*(durs-1000))); %tau*durs/max(durs);                  % Time constant of scaled data
maxProb = 0.5;

% Filter
filter_std = 25;                                % Standard deviation of the Gaussian filter used to estimate the underlying firing rate
Filter = @(t)(normpdf(t,0,filter_std));              % Spike filter
resolution = 1;
RemoveEdgeArtifact = 2*filter_std;

%% Generate activation functions
lambda = repmat(1-exp(-t/tau)',1,length(durs))*maxProb;       % Activation function of scaled data;
for i = 1:length(durs)
    lambdaScaled(:,i) = (1-exp(-t/taus(i))')*maxProb;
%    lambda(t > durs(i),i) = zeros(size(lambda(t>durs(i),i)));
%    lambdaScaled(t > durs(i),i) = zeros(size(lambdaScaled(t>durs(i),i)));
end

%% Generate spikes
% From coin flips
for i = 1:trials
    spikes(:,:,i) = rand(size(lambda)) < lambda;
    spikesScaled(:,:,i) = rand(size(lambdaScaled)) < lambdaScaled;
    for j = 1:length(durs)
        spikeT{1}{1,j}{i} = t(spikes(:,j,i)' & t <= durs(j));
        spikeTScaled{1}{1,j}{i} = t(spikesScaled(:,j,i)' & t <= durs(j));
    end
end

% Find rate by filtering spike train
for j = 1:length(durs)
    time = 0:resolution:durs(j);
    Spikes = horzcat(spikeT{1}{1,j}{:});
    SpikesScaled = horzcat(spikeTScaled{1}{1,j}{:});
    if isempty(Spikes)
        rates{1}{1}{j} = zeros(size(time));
    else
        Time = repmat(time',1,length(Spikes));
        Spikestemp = repmat(Spikes,size(Time,1),1);
        n = length(spikeT{1}{1,j});
        rates{1}{1}{j} = 1./n .* nansum( Filter(Time-Spikestemp) , 2);
    end
    if isempty(SpikesScaled)
        ratesDialated{1}{1}{j} = zeros(size(time));
    else
        Time = repmat(time',1,length(SpikesScaled));
        SpikesScaledtemp = repmat(SpikesScaled,size(Time,1),1);
        n = length(spikeTScaled{1}{1,j});
        ratesDialated{1}{1}{j} = 1./n .* nansum( Filter(Time-SpikesScaledtemp) , 2);
    end
end

% From ISI distribution p(ti+1 - ti) = lambda(ti+1)*exp(-?_(ti)^(ti+1)lambda(t')dt')
%spikeT = -repmat(lambda,[1 1 trials]).*log( rand(size(lambda,1),size(lambda,2),trials) );


%% Fit scaling factors to the data
[alpha, scaledTimes, scaledRates] = FitDialationFactor(spikeT{1}(1,:),durs,[1 1 1 1],'Filter',Filter,'resolution',resolution,'RemoveEdgeArtifact',RemoveEdgeArtifact);
ALPHA = vertcat(alpha);
ALPHA(ALPHA < 0) = NaN;

[alphaS, scaledTimesS, scaledRatesS] = FitDialationFactor(spikeTScaled{1}(1,:),durs,[1 1 1 1],'Filter',Filter,'resolution',resolution,'RemoveEdgeArtifact',RemoveEdgeArtifact);
ALPHAs = vertcat(alphaS);
ALPHAs(ALPHAs < 0) = NaN;