%% SimulateHMM
%
% Simulates a linear-Gaussian hidden Markov model with N dimensions and M
% measurements.
%
%%

%% Parameters
% General
N = 4;              % Dimensionality of the system
M = 50;            % Measurement dimensionality
T = 1000;           % Total number of time-stamps
trials = 50;       % Total number of trials
t = (1:T)';

% Transition model;
mu0 = [1 1 0 0];
Sig0 = [ 0.10 -0.10 0 0;...
        -0.10  0.10 0 0;...
         0     0    0 0;...
         0     0    0 0];
%A = [-0.15 -0.5; 0.5 -0.15];        % Deterministic portion
%SigT = [1 0; 0 1];      % Transition noise covariance
tau = 200;
SigT = [5 0 0 0;...
        0 5 0 0;...
        0 0 0.5 0;...
        0 0 0 0.5];
A = [1*tau  0 1  0;...
     0  1*tau 0  1;...
     -0.5*tau -0.7*tau 0  0;...
     0.7*tau  -0.5*tau 0  0];

% Measurement model
Beta = [randn(M,2) zeros(M,2)];%randn(M,N);
Gamma = zeros(M,N);
mu_c = -2;


% Plot options
nexamps = 5;    % Number of example trials to show in mean and variance plot
colors = [0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840;...
         0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840;...
         0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840;...
         0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840];
for i = 1:M
    if i <= size(colors,1)
        patchProperties{i}.FaceColor = colors(i,:) + 0.5*(ones(1,3)-colors(i,:));
    else
        break
    end
end

%% Simuate Hidden markov model
% Initialize the matrices
lambda = nan(T,M,trials);        % Measurements;
Y = nan(T,M,trials);        % Measurements;
X = nan(T,N,trials);        % Hidden state
dX = nan(T,N,trials);        % Change in hidden state

for k = 1:trials
    X(1,:,k) = mu0 + mvnrnd(zeros(1,N),Sig0);
    lambda(1,:,k) = exp( mu_c + Beta*X(1,:,k)' );
end


% Run simulation
for k = 1:trials
    for i = 2:T
        % Determine change in X
        Noise(i,:,k) = mvnrnd(zeros(N,1),SigT);
        
        % Update X
        X(i,:,k) = (A*X(i-1,:,k)' + Noise(i,:,k)')/tau;
        
        % Generate spiking responses
        lambda(i,:,k) = exp( mu_c + Beta*X(i,:,k)' );
        
    end
    for uniti = 1:M
        Y(:,uniti,k) = smooth(double(lambda(:,uniti,k) > rand(size(lambda(:,uniti,k)))),50);
    end
end

%% Find the mean velocity as a function of state
[Ybar, Yhat, dYhat, V, D, n2, mdY, q] = visualizerPCA(Y);

%% Calculate some statistics
% Mean and variance of x(i)
Xbar = mean(X,3);
Xvar = var(X,[],3);
 
% Ybar = mean(Y,3);
Yvar = var(Y,[],3);

%% Perform PCA *** LEGACY CODE ***
% Sigma = cov(Ybar);
% [V, D] = eig(Sigma);
% D = flipud(diag(D))/sum(diag(D));
% V = fliplr(V);
% 
% YhatBar = Ybar*V;
% 
% for i = 1:trials
%     Yhat(:,:,i) = Y(:,:,i)*V;
% end
% dYhat = diff(Yhat,1);
% 
% states = linspace(min(Yhat(:)),max(Yhat(:)),40);
% [STATES{1}, STATES{2}] = meshgrid(states);
% [INDX{1}, INDX{2}] = meshgrid(1:length(states));
% indx = [INDX{1}(:) INDX{2}(:)];
% [~, BIN] = histc(Yhat(1:end-1,:,:),states,3);
% bin1 = BIN(:,1,:);
% bin1 = bin1(:);
% bin2 = BIN(:,2,:);
% bin2 = bin2(:);
% for indxi = 1:size(indx,1)
%     accept = bin1 == indx(indxi,1) & bin2 == indx(indxi,2);
%     dytemp1 = dYhat(:,1,:);
%     dytemp1 = dytemp1(:);
%     dytemp2 = dYhat(:,2,:);
%     dytemp2 = dytemp2(:);
%     n2(indxi,1) = sum(accept);
%     mdY(indxi,1) = mean(dytemp1(accept));
%     mdY(indxi,2) = mean(dytemp2(accept));
% %     q(indxi,1) = mean( 1/2 *  sum(abs([dytemp1(accept) dytemp2(accept)]).^2) );
%     q(indxi,1) = mean( 1/2 *  sum(abs([mdY(indxi,1) mdY(indxi,2)]).^2) );
% end
% 
%% Plot the output

% Plot a random subset of trials and the mean/var
figure('Name','Mean and var of X')
trialinds2 = ceil(rand(nexamps,1)*trials);
for j = 1:2
    subplot(2,1,j)
    h = myPatch(t,Xbar(:,j),sqrt(Xvar(:,j)),'patchProperties',patchProperties{j});
    hold on
    plot(t,squeeze(X(:,j,trialinds2)),'Color',colors(j,:))
    %plot(t,squeeze(Y(:,1,trialinds2(end))),'o','Color',colors(j,:))
    
    plot(t,Xbar(:,j),'Color',colors(j,:),'LineWidth',3)
    xlabel('t')
    ylabel(['x_' num2str(j)])
    mymakeaxis(gca)
end

figure('Name','Example trajectories, first two dims')
plot(squeeze(X(:,1,trialinds2)),squeeze(X(:,2,trialinds2)),'.-','Color',[0.6 0.6 0.6]);
hold on
plot(squeeze(X(1,1,trialinds2)),squeeze(X(1,2,trialinds2)),'o','Color',[0.6 0.6 0.6]);
plot(squeeze(X(end,1,trialinds2)),squeeze(X(end,2,trialinds2)),'s','Color',[0.6 0.6 0.6]);
plot(Xbar(:,1),Xbar(:,2),'.-','Color',[0 0 0],'LineWidth',2,'MarkerSize',10);
plot(Xbar(1,1),Xbar(1,2),'o','Color',[0 0 0],'MarkerSize',10);
plot(Xbar(end,1),Xbar(end,2),'s','Color',[0 0 0],'MarkerSize',10);
axis equal
axis square
legend('E[x]','Example x')
xlabel('x_1')
ylabel('x_2')
mymakeaxis(gca)

figure('Name','Mean and var of Y')
trialinds2 = ceil(rand(nexamps,1)*trials);
for j = 1:3
    subplot(3,1,j)
    h = myPatch(t,Ybar(:,j),sqrt(Yvar(:,j)),'patchProperties',patchProperties{j});
    hold on
    
    plot(t,Ybar(:,j),'Color',colors(j,:),'LineWidth',3)
    axis tight
    xlabel('t')
    ylabel(['y_' num2str(j)])
    mymakeaxis(gca)
end

% PCA *** LEGACY CODE ***
% figure('Name','PCA')
% subplot(1,2,1)
% contour(STATES{1},STATES{2},reshape(n2,size(STATES{1})),20)
% colormap cool
% hold on
% allstates = [STATES{1}(:) STATES{2}(:)];
% cutoff = 0.001;
% quiver(allstates(n2/sum(n2) > cutoff,1),allstates(n2/sum(n2) > cutoff,2),...
%     mdY(n2/sum(n2) > cutoff,1),mdY(n2/sum(n2) > cutoff,2),2,'Color',[0 0 0])
% plot(YhatBar(:,1),YhatBar(:,2),'.-','LineWidth',1,'Color',[1 0 0])
% plot(YhatBar(1,1),YhatBar(1,2),'o','MarkerSize',10,'Color',[1 0 0])
% plot(YhatBar(end,1),YhatBar(end,2),'s','MarkerSize',10,'Color',[1 0 0])
% axis([-3 3 -3 3])
% axis square
% xlabel('PC_1')
% ylabel('PC_2')
% mymakeaxis(gca)
% 
% subplot(1,2,2)
% plot(cumsum(D),'ko')
% axis square
% xlabel('PCs')
% ylabel('Cumulative variance')
% mymakeaxis(gca)
% 
% figure('Name','q')
% surf(STATES{1},STATES{2},log(reshape(q,size(STATES{1}))),'EdgeColor','none')
% colormap cool
