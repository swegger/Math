%% BrodyTimerWrapper
%
%%

%% Variables

inNoiseLevels = linspace(0.00005,0.0002,10);

%% Run in loop
parfor i = 1:length(inNoiseLevels)
    [tp{i}, u{i}, v{i}, ncs{i}, thresholds{i}] = ...
        TimingWithBrodyNetwork('inNoise',inNoiseLevels(i));
end

%% Saving
save /om/user/swegger/analysis/Testing/BrodyWrapperTest