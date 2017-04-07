function [tp, u, v, ncs, threshold] = ...
    MagnitudeEstimatesWithBrodyNetwork(varargin)
%% TimingWithBrodyNetwork
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'dss',13:1:17)
addParameter(Parser,'N',3)
addParameter(Parser,'trialN',100)
addParameter(Parser,'t',0:400)
addParameter(Parser,'tau',80)
addParameter(Parser,'wI',0.0012)
addParameter(Parser,'uinit',1.5)
addParameter(Parser,'vinit',1.5)
addParameter(Parser,'gE',0.0015)
addParameter(Parser,'inNoise',0.0001)
addParameter(Parser,'transitionNoise',0.0001)
addParameter(Parser,'plotflg',false)

parse(Parser,varargin{:})

dss = Parser.Results.dss;
N = Parser.Results.N;
trialN = Parser.Results.trialN;
t = Parser.Results.t;
tau = Parser.Results.tau;
wI = Parser.Results.wI;
uinit = Parser.Results.uinit;
vinit = Parser.Results.vinit;
gE = Parser.Results.gE;
inNoise = Parser.Results.inNoise;
transitionNoise = Parser.Results.transitionNoise;
plotflg = Parser.Results.plotflg;

if length(uinit) == 1 && length(dss) ~=1
    uinit = uinit*ones(size(dss));
end
if length(vinit) == 1 && length(dss) ~=1
    vinit = vinit*ones(size(dss));
end
if length(gE) == 1 && length(dss) ~=1
    gE = gE*ones(size(dss));
end


tsColors = projectColorMaps('ts','samples',1:length(dss));

%% Iterate model for each d ss, no noise
for i = 1:N
    disp(['Epoch #' num2str(i)])
    if i == 1
        uic(i,:,:) = repmat(uinit,[1 1 trialN]);
        vic(i,:,:) = repmat(vinit,[1 1 trialN]);
        gEi(i,:,:) = repmat(mean(dss)/10000,[1 length(dss) trialN]);
    else
        uic(i,:,:) = uic(i-1,:,:) + ...
            (endpointsMeasure(i-1,:,:) - mean(endpoint(i-1,:)))/i;
        vic(i,:,:) = vic(i-1,:,:) + ...
            (endpointsMeasure(i-1,:,:) - mean(endpoint(i-1,:)))/i;
        gEi(i,:,:) = gEi(i-1,:,:) + ...
            (endpointsMeasure(i-1,:,:) - mean(endpoint(i-1,:)))/(i*1000);
    end
    
    for j = 1:length(dss)
        
        % Trying different intitial conditions
        [u{i,j}, v{i,j}, du{i,j}, dv{i,j}, ncs{i,j}] = BrodyNetwork(...
            squeeze(uic(i,j,:)),squeeze(vic(i,j,:)),'inputNoise',inNoise,...
            'trialN',trialN,'gE',squeeze(gEi(i,j,:)),...
            'transitionNoise',transitionNoise,...
            'plotflg',false,'t',t,'tau',tau,'wI',wI);
        
        [uMeasure{i,j}, vMeasure{i,j}, duMeasure{i,j}, dvMeasure{i,j}, ncsMeasure{i,j}] = ...
            BrodyNetwork(...
            squeeze(uic(i,j,:)),squeeze(vic(i,j,:)),'inputNoise',inNoise,...
            'trialN',trialN,'gE',dss(j)/10000,...
            'transitionNoise',transitionNoise,...
            'plotflg',false,'t',t,'tau',tau,'wI',wI);
        
    end
    
    % Calculate threshold
    temp = cat(3,v{i,:});
    tempMeasure = permute(cat(3,vMeasure{i,:}),[1 3 2]);
    endpoints = temp(end,:,:);
    endpointsMeasure(i,:,:) = tempMeasure(end,:,:);
    endpoint(i,:) = squeeze(mean(endpoints,2));
    threshold(i) = mean(endpoints(:));
    
end

 clear uic vic gEi

%% Run with input noise
for i = 1:N
    disp(['Epoch #' num2str(i)])
    if i == 1        
        uic(i,:,:) = repmat(uinit,[1 1 trialN]);
        vic(i,:,:) = repmat(vinit,[1 1 trialN]);
        gEi(i,:,:) = repmat(gE,[1 1 trialN]);
    else
        Eproj = errs(i-1,:,:);
        predErr(i-1,:) = mean(Eproj,3);
        uic(i,:,:) = repmat(uinit,[1 1 trialN]) - Eproj;
        vic(i,:,:) = repmat(vinit,[1 1 trialN]) - Eproj;
        gEi(i,:,:) = gEi(i-1,:,:) - Eproj/(i*1000);
    end
    
    for j = 1:length(dss)
        
        % Trying different intitial conditions
        [u{i,j}, v{i,j}, du{i,j}, dv{i,j}, ncs{i,j}] = BrodyNetwork(...
            squeeze(uic(i,j,:)),squeeze(vic(i,j,:)),'inputNoise',inNoise,...
            'trialN',trialN,'gE',squeeze(gEi(i,j,:)),...
            'transitionNoise',transitionNoise,...
            'plotflg',false,'t',t,'tau',tau,'wI',wI);
        
    end
    
    % Rotatate into "input" and "recurrent" dimensions, measure threshold
    % crossings and errors
    R = [cos(pi/4) -sin(pi/4); sin(pi/4) cos(pi/4)];
    for j = 1:length(dss)
        for triali = 1:trialN
            invec = [u{i,j}(:,triali),v{i,j}(:,triali)];
            outvec = R*invec';
            recurrentDim{i,j}(:,triali) = outvec(1,:);
            inputDim{i,j}(:,triali) = outvec(2,:);
            
        end
        
        % Calculate putative tps
        temp1 = repmat(t(:),[1,trialN]);
        temp2 = recurrentDim{i,j} > threshold(i);
        temp3 = temp1.*temp2;
        tp(i,j,:) = max(temp3,[],1);                % maximum t before all ts less than threshold for each trial
        
        % Calculate prediction errors
        errs(i,j,:) = recurrentDim{i,j}(t == dss(j),:) - threshold(i);
        
    end
    
end

Eproj = squeeze(errs(i,:,:));
predErr(N,:) = mean(Eproj,2);
% 
% for j = 1:length(tss)
%     us(j) = mean(u{end,j}(end,:));
%     vs(j) = mean(v{end,j}(end,:));
% end
% Eproj = [1 -1] * [us; vs];
% predErr(N,:) = (Eproj - mean(Eproj));

%% Evaluate dynamics on a grid
% us = linspace(0,4,20);
% vs = linspace(0,4,20);
% [Us,Vs] = meshgrid(us,vs);
% Us = Us(:);
% Vs = Vs(:);
% for j = 1:length(tss)
%     
%     % Trying different intitial conditions
%     [DU{j}, DV{j}] = BrodyStateSpaceDynamics(Us,Vs,...
%         'gE',gE(j),'plotflg',false,'t',0:tss(j),'tau',tau,'wI',wI);
%     
% end

%% Plotting
if plotflg
    
    % State space
    figure('Name','State space trajectories')
    for i = 1:N
        subplot(1,N,i)
        for j = 1:length(dss)
            plot(ncs{i,length(dss)-j+1}.x,ncs{i,length(dss)-j+1}.v,...
                'Color',tsColors(length(dss)-j+1,:) + 0.5*(ones(1,3)-tsColors(length(dss)-j+1,:)))
            hold on
            plot(ncs{i,length(dss)-j+1}.u,ncs{i,length(dss)-j+1}.x,...
                'Color',tsColors(length(dss)-j+1,:) + 0.5*(ones(1,3)-tsColors(length(dss)-j+1,:)))
            plot(1000*wI*u{i,length(dss)-j+1},...
                1000*wI*v{i,length(dss)-j+1},'.-',...
                'Color',tsColors(length(dss)-j+1,:))
        end
        recurr = threshold(i)*ones(1,2);
%         ins = vertcat(inputDim{i,:});
%         inp = [min(ins(:)) max(ins(:))];
%         Rback = R';%[cos(-pi/4) -sin(-pi/4); sin(-pi/4) cos(-pi/4)];
%         outs = Rback*[recurr; inp];
%         plot(outs(1,:),outs(2,:),'k--')
        axis square
        axis([0 4 0 4])
        xlabel('u')
        ylabel('v')
        title(['N = ' num2str(i)])
    end
    
%     figure('Name','State over time')
%     for j = 1:length(tss)
%         subplot(1,2,1)
%         plot(0:tss(length(tss)-j+1),u{length(tss)-j+1},...
%             'Color',tsColors(length(tss)-j+1,:))
%         hold on
%         
%         subplot(1,2,2)
%         plot(0:tss(length(tss)-j+1),v{length(tss)-j+1},...
%             'Color',tsColors(length(tss)-j+1,:))
%         hold on
%     end
%     subplot(1,2,1)
%     xlabel('t')
%     ylabel('u')
%     subplot(1,2,2)
%     xlabel('t')
%     ylabel('v')
%     
    % Prediction errors
    figure('Name','Prediction errors by number of samples')
    for i = 1:N
        plot(dss,predErr(i,:),'o-',...
            'Color',[0.8 0.8 0.8] - (i-1)*(0.8/(N-1)))
        hold on
        legendNames{i} = num2str(i);
    end
    legend(legendNames)
    plotHorizontal(0);
    xlabel('t_s')
    ylabel('Error')
    mymakeaxis(gca);
    
%     % State space dynamics
%     figure('Name','State space dynamics')
%     for j = 1:length(tss)
%         subplot(1,length(tss),j)
%         quiver(1000*wI*Us,1000*wI*Vs,DU{j},DV{j},...
%             'Color',tsColors(j,:) + 0.5*(ones(1,3)-tsColors(j,:)))
%         hold on
%         plot(ncs{j}.x,ncs{j}.v,...
%             'Color',tsColors(j,:) + 0.5*(ones(1,3)-tsColors(j,:)))
%         hold on
%         plot(ncs{j}.u,ncs{j}.x,...
%             'Color',tsColors(j,:) + 0.5*(ones(1,3)-tsColors(j,:)))
%         plot(1000*wI*u{j},1000*wI*v{j},'.-',...
%             'Color',tsColors(j,:))
%         
%         
%         axis square
%         axis([0 4 0 4])
%         xlabel('u')
%         ylabel('v')
%     end
    
    
end