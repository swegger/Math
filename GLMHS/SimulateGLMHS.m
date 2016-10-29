%% SimulateGLMHS
%
%   Simulates a hidden markov model with poisson measurements stemming from
%   an activation function that acts as a GLM.
%
%%

%% Parameters
% General
N = 2;              % Dimensionality of the system
M = 20;            % Measurement dimensionality
T = 1000;           % Total number of time-stamps
trials = 5;       % Total number of trials
t = (1:T)';

% Transition model
tau = 500;
A = [1 1;...
    -1/tau 0];
%A = [0 -1; 1 0];        % Deterministic portion
% A = [1  0     1  0;...
%      0  1     0  1;...
%      0 -1/tau 0  0;...
%      1/tau  0 0  0];

Gamma = [0.5/tau 0;...
         0 0.5/tau];
% Gamma = [1 0; 0 1];      % Transition noise covariance
% Gamma = [1  0  0  0;...
%          0  1  0  0;...
%          0  0  1  0;...
%          0  0  0  1]/tau;

% Measurement model
C = repmat([1 0],[M,1]);
%C = randn(M,N)/10;
%C(:,end-N/2+1:end) = 0;
mu_c = -1;


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
for i = 1:min([M size(colors,1)])
    patchProperties{i}.FaceColor = colors(i,:) + 0.5*(ones(1,3)-colors(i,:));
end
plotflg = false;

%% Simuate Hidden markov model
% Initialize the matrices
lambda = nan(M,T,trials);        % Measurements;
X = nan(M,T,trials);        % Measurements;
Z = nan(N,T,trials);        % Hidden state
%dZ = nan(T,N,trials);        % Change in hidden state

Z(:,1,:) = repmat(ones(size(Z(:,1,1))),[1 1 trials]);
Z(2,1,:) = 0;
for k = 1:trials
    lambda(:,1,k) = exp( mu_c + C*Z(:,1,k) );
end


%Z = repmat(ones(T,size(Z,2),1),[1 1 trials]);
% Run simulation
for k = 1:trials
    for i = 2:T
        Z(:,i,k) = A*Z(:,i-1,k) + mvnrnd([0;0],Gamma)'; %mGamma*randn(N,1);
%         % Determine change in Z
%         dZ(i,:,k) = permute( A*permute(Z(i-1,:,k),[2 3 1])...
%             + Gamma*randn(N,1), [3 1 2] )/tau;
%         
%         % Update Z
%         Z(i,:,k) = Z(i-1,:,k) + dZ(i,:,k);
        
        % Generate spiking responses
        lambda(:,i,k) = exp( mu_c + C*Z(:,i,k) );
        
    end
end
X = double(lambda > rand(size(lambda)));
%X = poissrnd(lambda);
% X(X > 1) = 1;

%% Fit model to the data
model0.A = A+randn(size(A))/10;
model0.Gamma = Gamma;%+randn(size(Gamma))/10;
% model0.Gamma = (model0.Gamma + model0.Gamma')/2;
model0.mu_c = mu_c*ones(M,1) +randn(M,1)/10; 
model0.C = C+randn(size(C))/10;
[model, LL, z, V, Vj, Vhat, V_, muhat, mu_] = ...
    identifyGLMHS_EM(X,model0,'mu0',Z(:,1,:),'verbose',true);

%% Calculate some statistics
% Mean and variance of x(i)
Zbar = mean(Z,3);
Zvar = var(Z,[],3);

Xbar = mean(X,3);%mean(Y,3);
Xvar = var(X,[],3);%var(Y,[],3);

%% Perform PCA
Sigma = cov(Xbar');
[Q, D] = eig(Sigma);
D = flipud(diag(D))/sum(diag(D));
Q = fliplr(Q);

YhatBar = Xbar'*Q;

for i = 1:trials
    Yhat(:,:,i) = X(:,:,i)'*Q;
end
dYhat = diff(Yhat,1);

states = linspace(min(Yhat(:)),max(Yhat(:)),80);
[STATES{1}, STATES{2}] = meshgrid(states);
[INDX{1}, INDX{2}] = meshgrid(1:length(states));
indx = [INDX{1}(:) INDX{2}(:)];
[~, BIN] = histc(Yhat(1:end-1,:,:),states,3);
bin1 = BIN(:,1,:);
bin1 = bin1(:);
bin2 = BIN(:,2,:);
bin2 = bin2(:);
for indxi = 1:size(indx,1)
    accept = bin1 == indx(indxi,1) & bin2 == indx(indxi,2);
    dytemp1 = dYhat(:,1,:);
    dytemp1 = dytemp1(:);
    dytemp2 = dYhat(:,2,:);
    dytemp2 = dytemp2(:);
    n2(indxi,1) = sum(accept);
    mdY(indxi,1) = mean(dytemp1(accept));
    mdY(indxi,2) = mean(dytemp2(accept));
end
mdY2 = mdY;
mdY2(n2/sum(n2) < 0.0001,:) = NaN;

%% Fit to the mean
% Theta0 = [1 0 reshape(A,1,numel(A)) reshape(Beta,1,numel(Beta))]...
%     + 0.000001*randn(size([1 0 reshape(A,1,numel(A)) reshape(Beta,1,numel(Beta))]));
% constraintA = [zeros(1,N) reshape(diag(ones(1,N)),1,numel(diag(ones(1,N)))) zeros(1,N*M)];
% constraintb = 0;        % Constrain trace(A) <= 0
% lowerBound = [Theta0(1:N)-0.01 -ones(1,N^2)/10 -inf(1,N*M)];
% upperBound = [Theta0(1:N)+0.01 ones(1,N^2)/10 inf(1,N*M)];
% [x0, Ahat, Betahat, Xhat, Yhat] = fitLDSmean(Ybar,N,'Theta0',Theta0,...
%     'constraintA',constraintA,'constraintb',constraintb,...
%     'lowerBound',lowerBound,'upperBound',upperBound);

%% Plot the output
if plotflg
    % Hidden States and measurements across trials
    figure('Name','Hidden states and measurements')
    trialinds = ceil(rand(2,1)*trials);
    for j = 1:N
        subplot(N,1,j)
        plot(squeeze(X(j,:,trialinds)),'o','Color',[0.7 0.7 0.7])
        hold on
        plot(squeeze(Z(j,:,trialinds)))
    end
    
    % Plot a random subset of trials and the mean/var
    figure('Name','Mean and var of X')
    trialinds2 = ceil(rand(nexamps,1)*trials);
    for j = 1:N
        subplot(N,1,j)
        h = myPatch(t,Zbar(j,:)',sqrt(Zvar(j,:))','patchProperties',patchProperties{j});
        hold on
        plot(t,squeeze(Z(j,:,trialinds2)),'Color',colors(j,:))
        %plot(t,squeeze(Y(:,1,trialinds2(end))),'o','Color',colors(j,:))
        
        plot(t,Zbar(j,:),'Color',colors(j,:),'LineWidth',3)
    end
    
    
    figure('Name','Mean and var of X')
    trialinds2 = ceil(rand(nexamps,1)*trials);
    for j = 1:3
        subplot(3,1,j)
        h = myPatch(t,Xbar(j,:)',sqrt(Xvar(j,:))','patchProperties',patchProperties{j});
        hold on
        %plot(t,squeeze(X(j,:,trialinds2)),'Color',colors(j,:))
        %plot(t,squeeze(Y(:,1,trialinds2(end))),'o','Color',colors(j,:))
        
        plot(t,Xbar(j,:),'Color',colors(j,:),'LineWidth',3)
    end
    
    figure('Name','Inferred state of Z from X, PCA')
    for j = 1:N
        subplot(N,1,j)
        plot(t,YhatBar(:,j)/max(YhatBar(:,j)),'LineWidth',3,'Color',[0.6 0.6 0.6])
        hold on
        plot(t,Zbar(j,:)/max(Zbar(j,:)),'LineWidth',3,'Color',[0 0 0])
    end
    
    figure('Name','Eigenvalues from PCA')
    plot(cumsum(D),'ko')
end

%% Save
save(['/om/user/swegger/analysis/SimulationResults/SimGLMHS_Test' datestr(now,'yyyymmdd')])
