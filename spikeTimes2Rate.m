function [rate, time, rateVariance, r] = spikeTimes2Rate(spikeTimes,varargin)
%% spikeTimes2Rate
%
%   [rate, time] = spikeTimes2Rate(spikeTimes);
%   Converts the spike times to an estimate of the underlying spiking
%   probability for each trial and a mean spiking probability across
%   trials.
%
%   [rate, time] = spikeTimes2Rate(spikeTimes,'Filter',filter)
%   Filters the spike trians according to the user supplied filter.
%
%   [rate...] = spikeTimes2Rate(spikeTimes,'StartEnd',StartEnd)
%   Defines a window to collect spikes and turn data into a spike train.
%
%   [rate...] = spikeTimes2Rate(...,'Histogram',HistogramOption)
%   Allows user to specify whether a histogram or filtering method is used.
%
%   [rate, time, rateVariance, r] =
%   spikeTimes2Rate(...,'ComputeVariance',ComputeVarianceOption)
%   Allow user to specify whether to compute the trial-by-trial rate (r) and
%   variance of the spike trains (rateVariance).
%
%
%%

% Defaults
defaultfilter = @(t)(double(~~t));

% Parse inputs
p = inputParser;
addRequired(p,'spikeTimes');
addParameter(p,'Filter',defaultfilter)
addParameter(p,'StartEnd',[]);
addParameter(p,'resolution',1/30);
addParameter(p,'Histogram','No')
addParameter(p,'ComputeVariance','No')

parse(p,spikeTimes,varargin{:})

spikeTimes = p.Results.spikeTimes;
Filter = p.Results.Filter;
StartEnd = p.Results.StartEnd;
resolution = p.Results.resolution;
Histogram = p.Results.Histogram;
ComputeVariance = p.Results.ComputeVariance;


% Collect times of all spikes
allspikes = [spikeTimes{:}];

% Determine the start and end times, if none supplied
if isempty(StartEnd)
    minT = min(allspikes);
    maxT = max(allspikes);
    time = [minT:resolution:maxT]';
else
    time = [StartEnd(1):resolution:StartEnd(2)]';
end

% Determine probability of firing
switch ComputeVariance
    case {'No','no','N','n',0}        
        % Find indices in time vector of each spike
        switch Histogram
            case {'No','no','N','n',0}
                % Find rate by filtering spike train
                if isempty(allspikes)
                    rate = zeros(size(time));
                else
                    Time = repmat(time,1,length(allspikes));
                    Spikes = repmat(allspikes,length(time),1);
                    rate = 1/length(spikeTimes) * sum( Filter(Time-Spikes) , 2);
                end
                r = nan(length(time),length(spikeTimes));
                rateVariance = nan(length(time),1);
                
            case {'Yes','yes','Y','y',1}
                % Find rate by binning across trials
                rate = histc(allspikes,time);
                rate = 1/length(spikeTimes) * rate(:);
                r = nan(length(time),length(spikeTimes));
                rateVariance = nan(length(time),1);
                
            otherwise
                error(['Histogram option ' Histogram ' not recognized!'])
        end
        
    case {'Yes','yes','Y','y',1}
        switch Histogram
            case {'No','no','N','n',0}
                % Find rate by filtering spike train
                for i = 1:length(spikeTimes)
                    if isempty(spikeTimes{i})
                        r(:,i) = zeros(size(time));
                    else
                        Time = repmat(time,1,length(spikeTimes{i}));
                        Spikes = repmat(spikeTimes{i},length(time),1);
                        r(:,i) = sum( Filter(Time-Spikes) , 2);
                    end
                end
                rate = mean(r,2);
                rateVariance = var(r,[],2);
            case {'Yes','yes','Y','y',1}
                % Find rate by binning across trials
                for i = 1:length(spikeTimes)
                    r(:,i) = histc(spikeTimes{i},time);
                end
                rate = mean(r,2);
                rateVariance = var(r,[],2);
                
            otherwise
                error(['Histogram option ' Histogram ' not recognized!'])
        end
    otherwise
        error(['Compute variance option ' ComputeVariance ' not recognized!'])
end