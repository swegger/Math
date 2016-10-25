function [ybar, yhat, dyhat, V, D, n, mdy, q] = visualizerPCA(y,varargin)
%% visualizerPCA
%
%   [ybar, yhat, dyhat, V, D, n, mdy, q] = visualizerPCA(y)
%
%       Performs PCA on input, y, and measures frequency of visits to each
%       state in the first 2 PCs and the mean change in state as a function
%       of state.
%       
%       y should be an TxMxN matrix where T is the number of time steps, M
%       is the number of neurons, and N is the number of trials.
%
%       Outputs the mean response (ybar), projections onto each PC (yhat),
%       change in response along each PC (dyhat), principal vectors (V),
%       eigenvalues (D), frequency of visits to each state (n), mean
%       velocity as a function of state (mdy), and a measure of the mean
%       speed as a function of state (q; see Sussillo and Barak, 2013).
%
%   ... = visualizerPCA(...,'states',states)
%       
%       Allows user to specify the states to measure. Defaults to 40
%       linearly spaced samples between the maximum and minimum projection 
%       values across all PCs.
%
%   ... = visualizerPCA(...,'nPCs',nPCs)
%
%       Allows user to specifiy the number of dimensions on which to
%       perform the analysis, using PCs 1 through nPCs as the dimensions.
%       Defaults to 2. ** NOTE ** Caution should be used here because
%       computational explosion with increasing nPCs.
%
%   ... = visualizerPCA(...,'cutoff',cutoff)
%
%       Allows user to specify the visualization of the velocity as a
%       function of the state by only showing velocities for those bins
%       where the frequency of visits was at least cutoff% of the data.
%       Defaults to 0.1% of the total number of data points. 
%
%
%   ... = visualizerPCA(...,'PlotOpts',PlotOpts)
%       
%       Allows user to control plotting.
%           PlotOpts.On - logical controling whether or not to plot;
%              defaults to true.
%           PlotOpts.PCsToPlot - 2 element vector allowing user to control
%              which PCs to plot. Defaults to [1 2] for the first 2 PCs.
%           PlotOpts.scale - controls scaling of velocity vectors in quiver
%              plot. Defaults to 2 time steps.
%
%
%   Version   by    date
%   1.0       swe   20161025
%%

%% Defaults
PlotOpts_default.On = true;
PlotOpts_default.PCsToPlot = [1 2];

%% Parse inputs

Parser = inputParser;

addRequired(Parser,'y')                                 % Input
addParameter(Parser,'states',NaN)                       % Vector of states to measure frequency of visits
addParameter(Parser,'nPCs',2)                           % Dimensionality of analysis
addParameter(Parser,'PlotOpts',PlotOpts_default)        % Controls plotting
addParameter(Parser,'cutoff',0.1)                    % Controls plotting of velocity states
addParameter(Parser,'scale',2)                          % Controls scaling of velocity vectors in plot

parse(Parser,y,varargin{:})

y = Parser.Results.y;
states = Parser.Results.states;
nPCs = Parser.Results.nPCs;
PlotOpts = Parser.Results.PlotOpts;
cutoff = Parser.Results.cutoff;

%% Perform PCA
ybar = mean(y,3);
Sigma = cov(ybar);
[V, D] = eig(Sigma);
D = flipud(diag(D))/sum(diag(D));
V = fliplr(V);

yhatBar = ybar*V;

yhat = nan(size(y));
for i = 1:size(y,3)
    yhat(:,:,i) = y(:,:,i)*V;
end
dyhat = diff(yhat,1);


%% Measure frequency of visits into each state and mean velocity as a function of state

% Set up list of states
if isnan(states)
    states = linspace(min(yhat(:)),max(yhat(:)),40);
end
STATEStemp = cell(1,nPCs);
[STATEStemp{:}] = ndgrid(states);
STATES = zeros(length(states)^nPCs,nPCs);
for j = 1:nPCs
    STATES(:,j) = [STATEStemp{j}(:)];
end
INDXtemp = cell(1,nPCs);
[INDXtemp{:}] = ndgrid(1:length(states));
INDX = zeros(length(states)^nPCs,nPCs);
for j = 1:nPCs
    INDX(:,j) = [INDXtemp{j}(:)];
end

% Label each time stamp with the bin number of that data
[~, BIN] = histc(yhat(1:end-1,:,:),states,3);
bin1 = BIN(:,1,:);
bin1 = bin1(:);
bin2 = BIN(:,2,:);
bin2 = bin2(:);
dytemp = reshape(permute(dyhat(:,1:nPCs,:),[1 3 2]),...
    [numel(dyhat(:,1:nPCs,:))/nPCs nPCs]);                  % Changes in state in PC space

% Cycle through state list to measure frequency of visits and mean velocity
l = size(INDX,1);
n = nan(l,1);
mdy = nan(l,nPCs);
q = nan(l,1);
for indxi = 1:l
    
    accept = bin1 == INDX(indxi,1) & bin2 == INDX(indxi,2);                 % visits to this state
    n(indxi,1) = sum(accept);                                              % Frequency of visits
    mdy(indxi,:) = mean(dytemp(accept,:),1);                                % mean velocity in state
    q(indxi,1) = mean( 1/2 *  sum(abs([mdy(indxi,1) mdy(indxi,2)]).^2) );   % mean speed   
    
end


%% Plot the output
if PlotOpts.On
        
    % Plot options
    if ~isfield(PlotOpts,'PCsToPlot')
        PlotOpts.PCsToPlot = [1 2];
    end
    if ~isfield(PlotOpts,'scale')
        PlotOpts.scale = 2;
    end
    
    % PCA
    figure('Name','PCA')
    subplot(1,2,1)
    contour(STATEStemp{PlotOpts.PCsToPlot(1)},STATEStemp{PlotOpts.PCsToPlot(2)},...
        reshape(n,size(STATEStemp{1})),20)
    colormap cool
    hold on
    
    quiver(STATES(n/sum(n) > cutoff/100,PlotOpts.PCsToPlot(1)),...
        STATES(n/sum(n) > cutoff/100,PlotOpts.PCsToPlot(2)),...
        mdy(n/sum(n) > cutoff/100,PlotOpts.PCsToPlot(1)),...
        mdy(n/sum(n) > cutoff/100,PlotOpts.PCsToPlot(2)),...
        PlotOpts.scale,'Color',[0 0 0])
    plot(yhatBar(:,1),yhatBar(:,2),'.-','LineWidth',1,'Color',[1 0 0])
    plot(yhatBar(1,1),yhatBar(1,2),'o','MarkerSize',10,'Color',[1 0 0])
    plot(yhatBar(end,1),yhatBar(end,2),'s','MarkerSize',10,'Color',[1 0 0])
    axis([-3 3 -3 3])
    axis square
    xlabel(['PC_' num2str(PlotOpts.PCsToPlot(1))])
    ylabel(['PC_' num2str(PlotOpts.PCsToPlot(2))])
    mymakeaxis(gca)
    
    subplot(1,2,2)
    plot(cumsum(D),'ko')
    axis square
    xlabel('PCs')
    ylabel('Cumulative variance')
    mymakeaxis(gca)
    
    figure('Name','q')
    surf(STATEStemp{PlotOpts.PCsToPlot(1)},STATEStemp{PlotOpts.PCsToPlot(2)},...
        log(reshape(q,size(STATEStemp{1}))),'EdgeColor','none')
    colormap cool
    xlabel(['PC_' num2str(PlotOpts.PCsToPlot(1))])
    ylabel(['PC_' num2str(PlotOpts.PCsToPlot(2))])
    zlabel('q')
    
end
