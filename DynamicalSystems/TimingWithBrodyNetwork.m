%% TimingWithBrodyNetwork
%
%
%%

%% Variables

% Task
tss = 600:100:1000;
N = 4;

% Network
tau = 240;
wI = 0.0012;
uinit = 1.5*ones(size(tss));
vinit = 1.5*ones(size(tss));
% uinit = fliplr(linspace(1.5,2.5,length(tss)));
% vinit = fliplr(linspace(1.5,2.5,length(tss)));
% uinit = linspace(1,2,length(tss));
% vinit = linspace(1,2,length(tss));

gE = 0.0015*ones(size(tss));
% gE = linspace(0.00125,0.00175,length(tss));

% Plotting
plotflg = false;
tsColors = projectColorMaps('ts','samples',1:length(tss));

%% Iterate model for each tss
for i = 1:N
    disp(['Epoch #' num2str(i)])
    if i == 1
        uic(i,:) = uinit;
        vic(i,:) = vinit;
        gEi(i,:) = gE;
    else
        for j = 1:length(tss)
            us(j) = u{i-1,j}(end);
            vs(j) = v{i-1,j}(end);
        end
        Eproj = [1 -1] * [us; vs];
        predErr(i-1,:) = (Eproj - mean(Eproj));
        uic(i,:) = uic(i-1,:) - (Eproj - mean(Eproj))/3;
        vic(i,:) = vic(i-1,:) - (Eproj - mean(Eproj))/3;
        gEi(i,:) = gEi(i-1,:) - (Eproj - mean(Eproj))/2000;
    end
    
    for j = 1:length(tss)
        
        % Trying different intitial conditions
        [u{i,j}, v{i,j}, du{i,j}, dv{i,j}, ncs{i,j}] = BrodyNetwork(...
            uic(i,j),vic(i,j),...
            'gE',gEi(i,j),'plotflg',false,'t',0:tss(j),'tau',tau,'wI',wI);
        
    end
end

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
            plot(ncs{i,length(tss)-j+1}.x,ncs{i,length(tss)-j+1}.v,...
                'Color',tsColors(length(tss)-j+1,:) + 0.5*(ones(1,3)-tsColors(length(tss)-j+1,:)))
            hold on
            plot(ncs{i,length(tss)-j+1}.u,ncs{i,length(tss)-j+1}.x,...
                'Color',tsColors(length(tss)-j+1,:) + 0.5*(ones(1,3)-tsColors(length(tss)-j+1,:)))
            plot(1000*wI*u{i,length(tss)-j+1},1000*wI*v{i,length(tss)-j+1},'.-',...
                'Color',tsColors(length(tss)-j+1,:))
        end
        axis square
        axis([0 4 0 4])
        xlabel('u')
        ylabel('v')
        title(['N = ' num2str(i)])
    end
    
    figure('Name','State over time')
    for j = 1:length(tss)
        subplot(1,2,1)
        plot(0:tss(length(tss)-j+1),u{length(tss)-j+1},...
            'Color',tsColors(length(tss)-j+1,:))
        hold on
        
        subplot(1,2,2)
        plot(0:tss(length(tss)-j+1),v{length(tss)-j+1},...
            'Color',tsColors(length(tss)-j+1,:))
        hold on
    end
    subplot(1,2,1)
    xlabel('t')
    ylabel('u')
    subplot(1,2,2)
    xlabel('t')
    ylabel('v')
    
    % Prediction errors
    figure('Name','Prediction errors by number of samples')
    for i = 1:N-1
        plot(tss,predErr(i,:),'o-',...
            'Color',[0.8 0.8 0.8] - (i-1)*(0.8/(N-2)))
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