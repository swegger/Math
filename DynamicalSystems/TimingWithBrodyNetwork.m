function [tp, u, v, ncs, threshold] = ...
    TimingWithBrodyNetwork(varargin)
%% TimingWithBrodyNetwork
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'tss',600:100:1000)
addParameter(Parser,'N',3)
addParameter(Parser,'trialN',1000)
addParameter(Parser,'t',0:1500)
addParameter(Parser,'tau',380)
addParameter(Parser,'wI',0.0012)
addParameter(Parser,'uinit',1.5)
addParameter(Parser,'vinit',1.5)
addParameter(Parser,'gE',0.0015)
addParameter(Parser,'inNoise',0.0001)
addParameter(Parser,'transitionNoise',0.0001)
addParameter(Parser,'plotflg',false)

parse(Parser,varargin{:})

tss = Parser.Results.tss;
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

if length(uinit) == 1 && length(tss) ~=1
    uinit = uinit*ones(size(tss));
end
if length(vinit) == 1 && length(tss) ~=1
    vinit = vinit*ones(size(tss));
end
if length(gE) == 1 && length(tss) ~=1
    gE = gE*ones(size(tss));
end

% Plotting
tsColors = projectColorMaps('ts','samples',1:length(tss));

%% Iterate model for each tss, no noise
for i = 1:N
    disp(['Epoch #' num2str(i)])
    if i == 1
        uic(i,:) = uinit;
        vic(i,:) = vinit;
        gEi(i,:) = gE;
    else
%         for j = 1:length(tss)
%             us(j) = mean(u{i-1,j}(t == tss(j),:));
%             vs(j) = mean(v{i-1,j}(t == tss(j),:));
%         end
%         Eproj = [1 -1] * [us; vs];
%         predErr(i-1,:) = (Eproj - mean(Eproj));
        Eproj = endpoint(i-1,:);
        uic(i,:) = uic(i-1,:) - (Eproj - mean(Eproj))/i;
        vic(i,:) = vic(i-1,:) - (Eproj - mean(Eproj))/i;
        gEi(i,:) = gEi(i-1,:) - (Eproj - mean(Eproj))/(i*1000);
    end
    
    for j = 1:length(tss)
        
        % Trying different intitial conditions
        [u{i,j}, v{i,j}, du{i,j}, dv{i,j}, ncs{i,j}] = BrodyNetwork(...
            squeeze(uic(i,j,:)),squeeze(vic(i,j,:)),'inputNoise',0,...
            'trialN',trialN,'gE',squeeze(gEi(i,j,:)),...
            'transitionNoise',0,...
            'plotflg',false,'t',t,'tau',tau,'wI',wI);
        
    end
    
    % Rotatate into "input" and "recurrent" dimensions, measure threshold
    % crossings and errors
    R = [cos(pi/4) -sin(pi/4); sin(pi/4) cos(pi/4)];
    for j = 1:length(tss)
        for triali = 1:trialN
            invec = [u{i,j}(:,triali),v{i,j}(:,triali)];
            outvec(:,:,triali) = R*invec';
        end
        endpoint(i,j) = mean(outvec(1,t == tss(j),:),3);
    end
    
    % Calculate threshold
    threshold(i) = mean(endpoint(i,:));
    
end

clear uic vic gEi outvec

%% Run with input noise
for i = 1:N
    disp(['Epoch #' num2str(i)])
    if i == 1        
        uic(i,:,:) = repmat(uinit,[1 1 trialN]);
        vic(i,:,:) = repmat(vinit,[1 1 trialN]);
        gEi(i,:,:) = repmat(gE,[1 1 trialN]);
    else
        Eproj = errs(i-1,:,:);
%         uic(i,:,:) = repmat(uinit,[1 1 trialN]) - Eproj;
%         vic(i,:,:) = repmat(vinit,[1 1 trialN]) - Eproj;
        uic(i,:,:) = uic(i-1,:,:) - Eproj/i;
        vic(i,:,:) = vic(i-1,:,:) - Eproj/i;
        gEi(i,:,:) = gEi(i-1,:,:) - Eproj/(i*1000);
    end
    
    for j = 1:length(tss)
        
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
    for j = 1:length(tss)
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
        errs(i,j,:) = recurrentDim{i,j}(t == tss(j),:) - threshold(i);
        
    end
    
end

Eproj = squeeze(errs(i,:,:));
predErr(N,:) = mean(Eproj,2);
predErr = squeeze(mean(errs,3));

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
        for j = 1:length(tss)
            samps = randsample(trialN,1);
            plot(ncs{i,length(tss)-j+1}.x(ncs{i,length(tss)-j+1}.x > 0),ncs{i,length(tss)-j+1}.v(samps,(ncs{i,length(tss)-j+1}.x > 0)),...
                'Color',tsColors(length(tss)-j+1,:) + 0.5*(ones(1,3)-tsColors(length(tss)-j+1,:)))
            hold on
            plot(ncs{i,length(tss)-j+1}.u(samps,(ncs{i,length(tss)-j+1}.x > 0)),ncs{i,length(tss)-j+1}.x(ncs{i,length(tss)-j+1}.x > 0),...
                'Color',tsColors(length(tss)-j+1,:) + 0.5*(ones(1,3)-tsColors(length(tss)-j+1,:)))
        end
        for j = 1:length(tss)
            plot(1000*wI*u{i,length(tss)-j+1}(t <= tss(length(tss)-j+1),samps),...
                1000*wI*v{i,length(tss)-j+1}(t <= tss(length(tss)-j+1),samps),'.-',...
                'Color',tsColors(length(tss)-j+1,:))
            hold on
        end
        recurr = threshold(i)*ones(1,2);
        ins = vertcat(inputDim{i,:});
        inp = [min(ins(:))-0.5 max(ins(:))+0.5];
        Rback = R';%[cos(-pi/4) -sin(-pi/4); sin(-pi/4) cos(-pi/4)];
        outs = Rback*[recurr; inp];
        plot(outs(1,:),outs(2,:),'k--')
        axis square
        axis([0 4 0 4])
        xlabel('u')
        ylabel('v')
        mymakeaxis(gca,'xytitle',['N = ' num2str(i)])
%         title(['N = ' num2str(i)])
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
        plot(tss,predErr(i,:),'o-',...
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