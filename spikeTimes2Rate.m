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
addParameter(p,'time',[]);

parse(p,spikeTimes,varargin{:})

spikeTimes = p.Results.spikeTimes;
Filter = p.Results.Filter;
StartEnd = p.Results.StartEnd;
resolution = p.Results.resolution;
Histogram = p.Results.Histogram;
ComputeVariance = p.Results.ComputeVariance;
time = p.Results.time;


% Collect times of all spikes
allspikes = [spikeTimes{:}];

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
                        for i = 1:length(spikeTimes)
                            Time = [Time repmat(time(:,i),1,length(spikeTimes{i}))];
                        end
                    end                        
                    Spikes = repmat(allspikes,size(Time,1),1);
                    n = sum(~isnan(time),2);
                    rate = 1./n .* nansum( Filter(Time-Spikes) , 2);
                end
                r = nan(size(time,1),length(spikeTimes));
                rateVariance = nan(size(time,1),1);
                
            case {'Yes','yes','Y','y',1}
                % Find rate by binning across trials
%                 if size(time,2) > 1
%                     error('Multiple time vector input for Histogram == Yes')
%                 end
                rate = histc(allspikes,time);
                rate = 1/length(spikeTimes) * rate(:);
                r = nan(size(time,1),length(spikeTimes));
                rateVariance = nan(size(time,1),1);
                
                % Calculate filter
                t0 = 0:resolution:max(max(time))/2;
                tvec = [-fliplr(t0(2:end)) t0];
                f = Filter(tvec);
                
                % Find trial by trial rates
                for i = 1:length(spikeTimes)
                    counttemp = histc(spikeTimes{i},time(~isnan(time(:,i)),i));
                    spikeCount(:,i) = [counttemp(:); zeros(sum(isnan(time(:,i))),1)];
                end
                
                % Calculate rate and variance
                n = sum(~isnan(time),2);
                rate = conv(sum(spikeCount,2)./n,f,'same');
                rateVariance = sum( (r - repmat(rate,1,size(r,2))).^2 .* 1./(repmat(n,1,size(r,2))-1) ,2);
                
            otherwise
                error(['Histogram option ' Histogram ' not recognized!'])
        end
        
    case {'Yes','yes','Y','y',1}
        switch Histogram
            case {'No','no','N','n',0}
                % Find rate by filtering spike train
                if size(time,2) == 1
                    time = repmat(time,1,length(spikeTimes));
                end
                for i = 1:length(spikeTimes)
                    if isempty(spikeTimes{i})
                        r(:,i) = zeros(size(time,1),1);
                    else
                        Time = repmat(time(:,i),1,length(spikeTimes{i}));
                        Spikes = repmat(spikeTimes{i},size(time,1),1);
                        r(:,i) = sum( Filter(Time-Spikes) , 2);
                         
%                         r(:,i) = zeros(size(time,1),1);
%                         for j = 1:length(spikeTimes{i})
%                             r(:,i) = r(:,i) + Filter(time(:,i) - spikeTimes{i}(j));
%                         end
                    end
                end
                % Calculate rate and variance
                n = sum(~isnan(time),2);
                rate = nansum(r.*1./repmat(n,1,size(r,2)),2);
                rateVariance = nansum( (r - repmat(rate,1,size(r,2))).^2 ,2) ./(n-1);
                %rate = mean(r,2);
                %rateVariance = var(r,[],2);
                
%                 THIS DOESN'T WORK RIGHT
%                 % Find all the spikes that aren't rejected by time matrix
%                 
%                 allspikes2 = [];
%                 for i = 1:length(spikeTimes)
%                     tcut = time(find(isnan(time(:,i)),1)-1,i);
%                     if ~isempty(tcut)
%                         allspikes2 = [allspikes2 spikeTimes{i}(spikeTimes{i} <= tcut)];
%                     else
%                         allspikes2 = [allspikes2 spikeTimes{i}];
%                     end
%                 end
%                 
%                 % Calculate the rate
%                 Time = min(time(:)):resolution:(size(time,1)-1)*resolution;
%                 for i = 1:length(Time)
%                     rate(i,1) = sum(Filter(Time(i)-allspikes2));
%                 end
%                 n = sum(~isnan(time),2);
%                 rate = rate./n;
%                 
%                 % Calculate the variance
%                 for i = 1:length(Time)
%                     rateVariance(i,1) = sum( (Filter(Time(i)-allspikes2) - rate(i)).^2 );
%                 end
%                 rateVariance = rateVariance./(n-1);
                
                
                
            case {'Yes','yes','Y','y',1}
                % Find rate by binning across trials
%                 if size(time,2) > 1
%                     error('Multiple time vector input for Histogram == Yes')
%                 end
                % Calculate filter
                t0 = 0:resolution:max(max(time))/2;
                tvec = [-fliplr(t0(2:end)) t0];
                f = Filter(tvec);
                
                % Find trial by trial rates
                for i = 1:length(spikeTimes)
                    rtemp = histc(spikeTimes{i},time(~isnan(time(:,i)),i));
                    r(:,i) = conv([rtemp(:); zeros(sum(isnan(time(:,i))),1)],f,'same');
                end
                
                % Calculate rate and variance
                n = sum(~isnan(time),2);
                rate = nansum(r.*1./repmat(n,1,size(r,2)),2);
                rateVariance = sum( (r - repmat(rate,1,size(r,2))).^2 .* 1./(repmat(n,1,size(r,2))-1) ,2);
                
            otherwise
                error(['Histogram option ' Histogram ' not recognized!'])
        end
    otherwise
        error(['Compute variance option ' ComputeVariance ' not recognized!'])
end