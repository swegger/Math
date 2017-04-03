%% BrodyTimerWrapper
%
%%

%% Variables
inNoiseLevels = linspace(0.00005,0.0002,10);
transNoiseLevels = linspace(0.5,1.5,10);
[cond1, cond2] = meshgrid(inNoiseLevels,transNoiseLevels);
NoiseLevels = [cond1(:) cond2(:)];
tss = 600:100:1000;
trialN = 1000;
N = 3;

plotflg = false;

%% Run in loop
parallelPool = parpool(20,'IdleTimeout',24*60);

parfor i = 1:size(NoiseLevels,1)
    [tp{i}, u{i}, v{i}, ncs{i}, thresholds{i}] = ...
        TimingWithBrodyNetwork('inNoise',NoiseLevels(i,1),...
        'transitionNoise',NoiseLevels(i,2));
end

%% Calculate BIAS/VARIANCE for each noise level
for i = 1:size(NoiseLevels,1)
    BIAS(i,:) = sqrt( mean((mean(tp{i},3) - repmat(tss,[N,1,1])).^2 ,2) );
    VAR(i,:) = mean( var(tp{i},[],3), 2);
end

%% Saving
save /om/user/swegger/analysis/Testing/BrodyWrapperTest

%% Clean up
delete(parallelPool)

%% Plotting
if plotflg
    plot(BIAS,sqrt(VAR),'o')
end