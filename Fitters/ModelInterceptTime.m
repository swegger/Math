function [responseT, invSp] = ModelInterceptTime( ...
    plotOn, ...
    weberParams, ...
    positions, ...
    nSimulations, ...
    invSpSupport, ...
    nInvSp)
%This is a program to simulate the tracking of a ball through a constant
%trajectory path, leading to an interception when the ball reaches a
%designated point. To simplify computations, and to keep the relationship
%between speed and timing distributions linear, I'm using inverse
%speed rather than speed. Also, we're assuming that there is no
%additional noise between the time interval measurement and the calculation
%of inverse speed. In this case, we can model the observer as directly
%taking measurements of inverse speed which are subject to the same
%scalar variance parameter (wm) as the time interval measurements. Thus it
%doesn't matter what the intercue distances are. As of now, I'm leaving
%them in the code though, in case I add an observer for which it matters.

%This model uses a uniform prior.

%If you want to use this to fit an RS+G task, just pass in a vector of ones
%for 'positions' of length n(sets) + 1.

if ~exist('plotOn','var')
    plotOn = 1;
end

%% Parameters
%Weber fractions for measurement and reproduction distributions
if ~exist('weberParams','var')
    weberMeasure = 0.2;
    weberProduce = 0.07;
else
    weberMeasure = weberParams(1);
    weberProduce = weberParams(2);
end

%Cue positions
if ~exist('positions','var')
    positions = [0 10 20 30];
    cueDist = diff(positions(1:end-1));
    targetDist = positions(end) - positions(end-1);
else
    cueDist = diff(positions(1:end-1));
    targetDist = positions(end) - positions(end-1);
end

%Number of simulations to run
if ~exist('nSimulations','var')
    nSimulations = 500;
end

%Domain, support, resolution, and actual values for inverse velocity.
%Domain is only used if you want to actually look at the posteriors.
if ~exist('invSpSupport','var')
    invSpSupport = [1/15 2/15];
end

if ~exist('nInvSp','var')
    nInvSp = 20;
end
dInvSp = 0.001;
domainLimsInvSp = [0.025 .2];
invSp = linspace(invSpSupport(1),invSpSupport(2),nInvSp);

%Alternate observer paramters. These were used as a figure in a Simons
%Postdoctoral Fellowship application, 2/28/2014. Only the standard ideal
%observer is an output of the simulation.
weberMeasure_NM = 1.68*weberMeasure;

weberMeasure_HP = 1.28*weberMeasure;
weberMeasureAssumed_HP = .15*weberMeasure_HP; %too low and the sim breaks.

weberProduce_MN = 1.75*weberProduce;

%% Internal Prior for inverse velocity
[fInvSp, domainInvSp] = ...
get_uniform_distribution(dInvSp,domainLimsInvSp,invSpSupport);

%% Estimates and errors
%Optimal observer
% measurementsT = zeros(nSimulations,numel(invSp),numel(cueDist));
calculatedInvSp = zeros(nSimulations,numel(invSp),numel(cueDist));
estInvSp = zeros(nSimulations,numel(invSp),numel(cueDist));
varEstInvSp = zeros(nSimulations,numel(invSp),numel(cueDist));
biasEstInvSp = zeros(nSimulations,numel(invSp),numel(cueDist));
mseEstInvSp = zeros(nSimulations,numel(invSp),numel(cueDist));

tTarget = zeros(nSimulations,numel(invSp));
responseT = zeros(nSimulations,numel(invSp));
varResponseT = zeros(numel(invSp),1);
biasResponseT = zeros(numel(invSp),1);
mseResponseT = zeros(numel(invSp),1);

%Noisy measurements (NM) but otherwise optimal
% measurementsT_NM = zeros(nSimulations,numel(invSp),numel(cueDist));
calculatedInvSp_NM = zeros(nSimulations,numel(invSp),numel(cueDist));
estInvSp_NM = zeros(nSimulations,numel(invSp),numel(cueDist));
varEstInvSp_NM = zeros(nSimulations,numel(invSp),numel(cueDist));
biasEstInvSp_NM = zeros(nSimulations,numel(invSp),numel(cueDist));
mseEstInvSp_NM = zeros(nSimulations,numel(invSp),numel(cueDist));

tTarget_NM = zeros(nSimulations,numel(invSp));
responseT_NM = zeros(nSimulations,numel(invSp));
varResponseT_NM = zeros(numel(invSp),1);
biasResponseT_NM = zeros(numel(invSp),1);
mseResponseT_NM = zeros(numel(invSp),1);

%Overweighting prediction or "hypoprior" (HP)
% measurementsT_HP = zeros(nSimulations,numel(invSp),numel(cueDist));
calculatedInvSp_HP = zeros(nSimulations,numel(invSp),numel(cueDist));
estInvSp_HP = zeros(nSimulations,numel(invSp),numel(cueDist));
varEstInvSp_HP = zeros(nSimulations,numel(invSp),numel(cueDist));
biasEstInvSp_HP = zeros(nSimulations,numel(invSp),numel(cueDist));
mseEstInvSp_HP = zeros(nSimulations,numel(invSp),numel(cueDist));

tTarget_HP = zeros(nSimulations,numel(invSp));
responseT_HP = zeros(nSimulations,numel(invSp));
varResponseT_HP = zeros(numel(invSp),1);
biasResponseT_HP = zeros(numel(invSp),1);
mseResponseT_HP = zeros(numel(invSp),1);

%High motor noise (MN) but other wise optimal
responseT_MN = zeros(nSimulations,numel(invSp));
varResponseT_MN = zeros(numel(invSp),1);
biasResponseT_MN = zeros(numel(invSp),1);
mseResponseT_MN = zeros(numel(invSp),1);


%% Run simulation
for i = 1:numel(invSp)
    for j = 1:numel(cueDist)
        %% Optimal observer
        calculatedInvSp(:,i,j) = ...
            randn(1,nSimulations)*weberMeasure*invSp(i) + invSp(i);
        
        %Get estimates and posterior
        estInvSp(:,i,j) = bayesian_estimate( ...
            domainInvSp, ...
            squeeze(calculatedInvSp(:,i,1:j)), ...
            weberMeasure, ...
            invSpSupport);
        
        %Variance and mean squared error
        varEstInvSp(i,j) = var(estInvSp(:,i,j));
        biasEstInvSp(i,j) = mean(estInvSp(:,i,j)) - invSp(i);
        mseEstInvSp(i,j) = mean((estInvSp(:,i,j) - invSp(i)).^2);
        
        %% Noisy measurements (NM) but otherwise optimal
        calculatedInvSp_NM(:,i,j) = ...
            randn(1,nSimulations)*weberMeasure_NM*invSp(i) + invSp(i);
        
        %Get estimates and posterior
        estInvSp_NM(:,i,j) = bayesian_estimate( ...
            domainInvSp, ...
            squeeze(calculatedInvSp_NM(:,i,1:j)), ...
            weberMeasure_NM, ...
            invSpSupport);
        
        %Variance and mean squared error
        varEstInvSp_NM(i,j) = var(estInvSp_NM(:,i,j));
        biasEstInvSp_NM(i,j) = mean(estInvSp_NM(:,i,j)) - invSp(i);
        mseEstInvSp_NM(i,j) = mean((estInvSp_NM(:,i,j) - invSp(i)).^2);
        
        %% Overweighting prediction or "hypoprior" (HP)
        calculatedInvSp_HP(:,i,j) = ...
            randn(1,nSimulations)*weberMeasure_HP*invSp(i) + invSp(i);
        
        %Get estimates and posterior
        estInvSp_HP(:,i,j) = bayesian_estimate( ...
            domainInvSp, ...
            squeeze(calculatedInvSp_HP(:,i,1:j)), ...
            weberMeasureAssumed_HP, ...
            invSpSupport);
        
        %Variance and mean squared error
        varEstInvSp_HP(i,j) = var(estInvSp_HP(:,i,j));
        biasEstInvSp_HP(i,j) = mean(estInvSp_HP(:,i,j)) - invSp(i);
        mseEstInvSp_HP(i,j) = mean((estInvSp_HP(:,i,j) - invSp(i)).^2);
    end
%% Response interval
tActual = invSp(i)*targetDist;
%Ideal observer
tTarget(:,i) = targetDist*estInvSp(:,i,end);
responseT(:,i) = ...
    randn(1,nSimulations)'.*tTarget(:,i)*weberProduce + tTarget(:,i);

varResponseT(i) = var(responseT(:,i));
biasResponseT(i) = mean(responseT(:,i)) - tActual;
mseResponseT(i) = mean((responseT(:,i) - tActual).^2);

%Noisy measurements (NM) but otherwise optimal
tTarget_NM(:,i) = targetDist*estInvSp_NM(:,i,end);
responseT_NM(:,i) = ...
    randn(1,nSimulations)'.*tTarget_NM(:,i)*weberProduce + tTarget_NM(:,i);

varResponseT_NM(i) = var(responseT_NM(:,i));
biasResponseT_NM(i) = mean(responseT_NM(:,i)) - tActual;
mseResponseT_NM(i) = mean((responseT_NM(:,i) - tActual).^2);

%Overweighting prediction or "hypoprior" (HP)
tTarget_HP(:,i) = targetDist*estInvSp_HP(:,i,end);
responseT_HP(:,i) = ...
    randn(1,nSimulations)'.*tTarget_HP(:,i)*weberProduce + tTarget_HP(:,i);

varResponseT_HP(i) = var(responseT_HP(:,i));
biasResponseT_HP(i) = mean(responseT_HP(:,i)) - tActual;
mseResponseT_HP(i) = mean((responseT_HP(:,i) - tActual).^2);

%High motor noise (MN) but other wise optimal
responseT_MN(:,i) = ...
    randn(1,nSimulations)'.*tTarget(:,i)*weberProduce_MN + tTarget(:,i);

varResponseT_MN(i) = var(responseT_MN(:,i));
biasResponseT_MN(i) = mean(responseT_MN(:,i)) - tActual;
mseResponseT_MN(i) = mean((responseT_MN(:,i) - tActual).^2);
end

if plotOn
    %% Plot time estimates for different types of observers
    % figure('Position',[680   849   628   249]);
    % subplot(1,2,1)
    figure('Position',[680   858   253   240])
    % Overweighting prediction or "hypoprior" (HP)
    patch([invSp invSp(end:-1:1)], ...
        [mean(responseT_HP) + varResponseT_HP' ...
        mean(responseT_HP(:,end:-1:1)) - varResponseT_HP(end:-1:1)'], ...
        [.8 .8 1])
    hold on;
    h_HP = plot(invSp,mean(responseT_HP), ...
        'color',[.3 .3 .8], ...
        'linewidth',2);
    
    % Noisy measurements (NM) but otherwise optimal
    patch([invSp invSp(end:-1:1)], ...
        [mean(responseT_NM) + sqrt(varResponseT_NM') ...
        mean(responseT_NM(:,end:-1:1)) - sqrt(varResponseT_NM(end:-1:1)')], ...
        [1 .8 .8])
    hold on;
    h_NM = plot(invSp,mean(responseT_NM), ...
        'color',[.8 .3 .3], ...
        'linewidth',2);
    
    % High motor noise (MN) but other wise optimal
    patch([invSp invSp(end:-1:1)], ...
        [mean(responseT_MN) + sqrt(varResponseT_MN') ...
        mean(responseT_MN(:,end:-1:1)) - sqrt(varResponseT_MN(end:-1:1)')], ...
        [.8 1 .8])
    hold on;
    h_MN = plot(invSp,mean(responseT_MN), ...
        'color',[.3 .8 .3], ...
        'linewidth',2);
    
    %Ideal observer
    patch([invSp invSp(end:-1:1)], ...
        [mean(responseT) + sqrt(varResponseT)' ...
        mean(responseT(:,end:-1:1)) - sqrt(varResponseT(end:-1:1))'], ...
        [.8 .8 .8])
    hold on;
    h = plot(invSp,mean(responseT), ...
        'color',[.3 .3 .3], ...
        'linewidth',2);
    
    %Correct times
    plot(invSp,invSp*targetDist, ...
        'color','k', ...
        'linewidth',2, ...
        'linestyle','--')
    
    legend([h,h_NM,h_HP,h_MN],{num2str(sqrt(mean(mseResponseT))), ...
        num2str(sqrt(mean(mseResponseT_NM))), ...
        num2str(sqrt(mean(mseResponseT_HP))), ...
        num2str(sqrt(mean(mseResponseT_MN)))}, ...
        'Location','Northwest')
    
    axis tight square
    % Plot Bias vs. SD for different observers
    % subplot(1,2,2)
    figure('Position',[680   858   253   240]);
    plot(sqrt(mean(biasResponseT.^2)),sqrt(mean(varResponseT)), ...
        'Marker','.', ...
        'MarkerSize',30, ...
        'Color',[.3 .3 .3])
    hold on;
    h_NM = plot(sqrt(mean(biasResponseT_NM.^2)),sqrt(mean(varResponseT_NM)), ...
        'Marker','.', ...
        'MarkerSize',30, ...
        'Color',[.8 .3 .3]);
    h_HP = plot(sqrt(mean(biasResponseT_HP.^2)),sqrt(mean(varResponseT_HP)), ...
        'Marker','.', ...
        'MarkerSize',30, ...
        'Color',[.3 .3 .8]);
    h_MN = plot(sqrt(mean(biasResponseT_MN.^2)),sqrt(mean(varResponseT_MN)), ...
        'Marker','.', ...
        'MarkerSize',30, ...
        'Color',[.3 .8 .3]);
    
    set(gca, 'XLim', [0 .25], ...
        'YLim', [0 .25]);
    
    % small RMSE arc
    angles = linspace(0,pi/2,15);
    bias = sqrt(mean(mseResponseT))*cos(angles);
    sqrtVariance = sqrt(mean(mseResponseT))*sin(angles);
    h_MSE = plot(bias,sqrtVariance, ...
        'Color', [.5 .5 .5], ...
        'LineStyle', '--', ...
        'LineWidth', 1.5);
    
    % increased RMSE arc
    meanRmsePoor = mean([sqrt(mean(mseResponseT_NM)) sqrt(mean(mseResponseT_HP)) ...
        sqrt(mean(mseResponseT_MN))]);
    angles = linspace(0,pi/2,15);
    bias = meanRmsePoor*cos(angles);
    sqrtVariance = meanRmsePoor*sin(angles);
    plot(bias,sqrtVariance, ...
        'Color', [.5 .5 .5], ...
        'LineStyle', '--', ...
        'LineWidth', 1.5)
    
    legend([h_NM,h_MN,h_HP],{'Measurement noise', ...
        'Motor noise', ...
        'Error over-weighting'})
    text(1.05*bias(5),1.05*sqrtVariance(5),'Error contours')
    
    
    set(gca,'Box','Off', ...
        'XTickLabel',[], ...
        'YTickLabel',[])
    xlabel('Bias')
    ylabel('\surdVar')
    
    axis square
end

%% Other Things to plot
% % %
% % figure;
% % subplot(3,1,1)
% % hist(calculatedInvSp)
% % set(gca,'XLim',domainLimsInvSp)
% % subplot(3,1,2)
% % plot(domainInvSp,priorInvSp)
% % set(gca,'XLim',domainLimsInvSp)



% %%
% figure;
% subplot(3,1,1)
% hist(responseT(:,1))
% subplot(3,1,2)
% hist(estInvSp(:,i,end))
% subplot(3,1,3)
% plot(invSp,mean(responseT))
% 
% %% Plot velocity estimates through cues
% figure;
% plot(invSp,squeeze(mean(estInvSp)), ...
%     'color',[.3 .3 .3], ...
%     'linewidth',2)



