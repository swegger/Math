function [rate, time, rateVariance, r, mspikeCount, spikeCount] = spikeTimes2Rate(spikeTimes,varargin)
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
addParameter(p,'time',[]);
addParameter(p,'mixedUnits',false)
addParameter(p,'units',NaN)

parse(p,spikeTimes,varargin{:})

spikeTimes = p.Results.spikeTimes;
Filter = p.Results.Filter;
StartEnd = p.Results.StartEnd;
resolution = p.Results.resolution;
Histogram = p.Results.Histogram;
ComputeVariance = p.Results.ComputeVariance;
time = p.Results.time;
mixedUnits = p.Results.mixedUnits;
units = p.Results.units;

% Collect times of all spikes
if mixedUnits
    allspikes = double(spikeTimes{1});
else
    allspikes = double([spikeTimes{:}]);
end

% Determine the start and end times, if none supplied
if isempty(time)
    if isempty(StartEnd)
        minT = min(allspikes);
        maxT = max(allspikes);
        time = (minT:resolution:maxT)';
    else
        time = (StartEnd(1):resolution:StartEnd(2))';
    end
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
                    if size(time,2) == 1
                        Time = repmat(time,1,length(allspikes));
                    else
                        Time = [];
                        if mixedUnits
                            if isnan(units)
                                units = unique(spikeTimes{2});
                            end
                            for uniti = 1:length(units)
                                Time = [Time repmat(time(:,uniti),1,sum(spikeTimes{2} == units(uniti)))];
                            end
                        else
                            for i = 1:length(spikeTimes)
                                Time = [Time repmat(time(:,i),1,length(spikeTimes{i}))];
                            end
                        end
                    end                        
                    Spikes = repmat(allspikes,size(Time,1),1);
                    n = sum(~isnan(time),2);
                    rate = 1./n .* nansum( Filter(Time-Spikes) , 2);
                end
                if mixedUnits
                    if isnan(units)
                        units = unique(spikeTimes{2});
                    end
                    r = nan(size(time,1),length(units));
                else
                    r = nan(size(time,1),length(spikeTimes));
                end
                rateVariance = nan(size(time,1),1);
                mspikecount = nan(size(time,1),1);
                spikeCount = nan(size(time));
                
                
            case {'Yes','yes','Y','y',1}                
                % Calculate filter
                t0 = 0:resolution:max(max(time))/2;
                tvec = [-fliplr(t0(2:end)) t0];
                f = Filter(tvec);
                
                % Find trial by trial rates
                if size(time,2) == 1
                    if mixedUnits
                        if isnan(units)
                            units = unique(spikeTimes{2});
                        end
                        time = repmat(time,1,length(units));
                    else
                        time = repmat(time,1,length(spikeTimes));
                    end
                end
                if mixedUnits
                    if isnan(units)
                        units = unique(spikeTimes{2});
                    end
                    for uniti = 1:length(units)
                        spks = spikeTimes{1}(spikeTimes{2} == units(uniti));
                        if isempty(spks)
                            spikeCount(:,i) = zeros(size(time(:,i)));
                        else
                            counttemp = histc(spks,time(~isnan(time(:,i)),1));
                            spikeCount(:,i) = [counttemp(:); zeros(sum(isnan(time(:,i))),1)];
                        end
                    end
                else
                    for i = 1:length(spikeTimes)
                        if isempty(spikeTimes{i})
                            spikeCount(:,i) = zeros(size(time(:,i)));
                        else
                            counttemp = histc(spikeTimes{i},time(~isnan(time(:,i)),i));
                            spikeCount(:,i) = [counttemp(:); zeros(sum(isnan(time(:,i))),1)];
                        end
                    end
                end
                
                % Calculate rate and variance
                n = sum(~isnan(time),2);
                mspikeCount = sum(spikeCount,2)./n;
                mspikeCount(isnan(mspikeCount)) = 0;
                rate = conv(mspikeCount,f,'same');
                if mixedUnits
                    if isnan(units)
                        units = unique(spikeTimes{2});
                    end
                    r = nan(size(time,1),length(units));
                else
                    r = nan(size(time,1),length(spikeTimes));
                end
                rateVariance = nan(size(time,1),1);
                
            otherwise
                error(['Histogram option ' Histogram ' not recognized!'])
        end
        
    case {'Yes','yes','Y','y',1}
        switch Histogram
            case {'No','no','N','n',0}
                % Find rate by filtering spike train
                if size(time,2) == 1
                    if mixedUnits
                        if isnan(units)
                            units = unique(spikeTimes{2});
                        end
                        time = repmat(time,1,length(units));
                    else
                        time = repmat(time,1,length(spikeTimes));
                    end
                end
                if mixedUnits
                    if isnan(units)
                        units = unique(spikeTimes{2});
                    end
                    for uniti = 1:length(units)
                        spks = spikeTimes{1}(spikeTimes{2} == units(uniti));
                        if isempty(spks)
                            r(:,uniti) = zeros(size(time,1),1);
                        else
                            Time = repmat(time(:,uniti),1,length(spks));
                            Spikes = repmat(spks,size(time,1),1);
                            r(:,uniti) = sum( Filter(Time-Spikes) , 2);
                            
                        end
                    end
                else
                    for i = 1:length(spikeTimes)
                        if isempty(spikeTimes{i})
                            r(:,i) = zeros(size(time,1),1);
                        else
                            Time = repmat(time(:,i),1,length(spikeTimes{i}));
                            Spikes = repmat(spikeTimes{i},size(time,1),1);
                            r(:,i) = sum( Filter(Time-Spikes) , 2);
                            
                        end
                    end
                end
                % Calculate rate and variance
                n = sum(~isnan(time),2);
                rate = nansum(r.*1./repmat(n,1,size(r,2)),2);
                rateVariance = nansum( (r - repmat(rate,1,size(r,2))).^2 ,2) ./(n-1);
                mspikecount = nan(size(time,1),1);
                spikeCount = nan(size(time));               
                
            case {'Yes','yes','Y','y',1}
                % Calculate filter
                t0 = 1:resolution:floor(size(time,1)/2);
                tvec = [-fliplr(t0(1:end-2)) 0 t0(1:end-2)];
                f = Filter(tvec);
                f = f(f > 10^-10);
                
                % Find trial by trial rates and counts
                if size(time,2) == 1
                    if mixedUnits
                        if isnan(units)
                            units = unique(spikeTimes{2});
                        end
                        time = repmat(time,1,length(units));
                    else
                        time = repmat(time,1,length(spikeTimes));
                    end
                end
                
                if mixedUnits
                    if isnan(units)
                        units = unique(spikeTimes{2});
                    end
                    for i = 1:length(units)
                        spks = spikeTimes{1}(spikeTimes{2} == units(i));
                        % Counts
                        counttemp = histc(spks,time(~isnan(time(:,i)),i));
                        if isempty(counttemp)
                            spikeCount(:,i) = zeros(size(time(:,i)));
                        else
                            spikeCount(:,i) = [counttemp(:); zeros(sum(isnan(time(:,i))),1)];
                        end
                        
                        % Rates
                        rtemp = histc(spks,time(~isnan(time(:,i)),i));
                        if isempty(rtemp)
                            r(:,i) = zeros(size(time,1),1);
                        else
                            r(:,i) = [conv(rtemp(:),f,'same'); nan(sum(isnan(time(:,i))),1)];
                        end
                    end                    
                else
                    for i = 1:length(spikeTimes)
                        % Counts
                        counttemp = histc(spikeTimes{i},time(~isnan(time(:,i)),i));
                        if isempty(counttemp)
                            spikeCount(:,i) = zeros(size(time(:,i)));
                        else
                            spikeCount(:,i) = [counttemp(:); zeros(sum(isnan(time(:,i))),1)];
                        end
                        
                        % Rates
                        rtemp = histc(spikeTimes{i},time(~isnan(time(:,i)),i));
                        if isempty(rtemp)
                            r(:,i) = zeros(size(time,1),1);
                        else
                            r(:,i) = [conv(rtemp(:),f,'same'); nan(sum(isnan(time(:,i))),1)];
                        end
                    end
                end
                
                % Calculate rate and variance
                n = sum(~isnan(time),2);
                mspikeCount = sum(spikeCount,2)./n;
                rate = conv(sum(spikeCount,2)./n,f,'same');         % Calculate rate from counts; note mean(r) and this can potentially give different outputs, if time vector has NaNs
                rateVariance = sum( (r - repmat(rate,1,size(r,2))).^2 .* 1./(repmat(n,1,size(r,2))-1) ,2);
                
            otherwise
                error(['Histogram option ' Histogram ' not recognized!'])
        end
    otherwise
        error(['Compute variance option ' ComputeVariance ' not recognized!'])
end