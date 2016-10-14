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
Sig0 = 0.25;
%A = [-0.15 -0.5; 0.5 -0.15];        % Deterministic portion
%SigT = [1 0; 0 1];      % Transition noise covariance
tau = 100;
SigT = [5 0 0 0;...
        0 5 0 0;...
        0 0 0.5 0;...
        0 0 0 0.5];
A = [0.999*tau  0 1  0;...
     0  0.999*tau 0  1;...
     0 -0.5*tau 0  0;...
     0.5*tau  0 0  0];

% For fitting
BetaSize = 100;
omega = linspace(0,2*pi,BetaSize);
BasisN = 10;
psi = BasisN/4*ones(1,BasisN);
phi = (0:pi/2:BasisN*pi/2);%linspace(0,2*pi,BasisN);
for i = 1:BasisN
    BasisSet.k(i,:) = zeros(1,length(omega));
    q = omega > (-pi/2 + phi(i))/psi(i) & omega < (pi/2 + phi(i))/psi(i);
    BasisSet.k(i,q) = ( cos( omega(q)*psi(i) - phi(i) ) );
end
BasisSet.h = BasisSet.k;

% Measurement model
Beta = [randn(M,2) zeros(M,2)];%randn(M,N);
Gamma = zeros(M,N);
mu_c = -2;
% Beta = 0.01*flipud(BasisSet.k(2,:)');
% Gamma = -0.005*flipud(BasisSet.h(1,:)');
% GammaSize = length(Gamma);
% mu_c = 0;
% BetaParams = rand(N,M)/2;
% Beta = zeros(BetaSize,N,M);
% for i = 1:N
%     for j = 1:M
%         Beta(:,i,j) = exp(-BetaParams(i,j)*(0:BetaSize-1))/sum(exp(-BetaParams(i,j)*(0:BetaSize-1)));
%     end
% end
% Beta = flipud(Beta);
% %Beta = reshape(Beta,[BetaSize*N, M]);
% GammaSize = 100;
% GammaParams = rand(M,1);
% Gamma = repmat((0:GammaSize-1)',[1, M]);
% for j = 1:M
%     Gamma(:,j) = -exp(-GammaParams(j)*Gamma(:,j))/sum(exp(-GammaParams(j)*Gamma(:,j)));
% end
% Gamma = flipud(Gamma);

% Beta = 1;
% Gamma = -0.2;


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
    patchProperties{i}.FaceColor = colors(i,:) + 0.5*(ones(1,3)-colors(i,:));
end

%% Simuate Hidden markov model
% Initialize the matrices
lambda = nan(T,M,trials);        % Measurements;
Y = nan(T,M,trials);        % Measurements;
X = nan(T,N,trials);        % Hidden state
dX = nan(T,N,trials);        % Change in hidden state

X(1,:,:) = 1*repmat(ones(size(X(1,:,1))),[1 1 trials]) + randn(size(X(1,:,:)))*Sig0;
% tempBeta = reshape(Beta(end,:,:),[N, M])';
for k = 1:trials
    lambda(1,:,k) = exp( mu_c + Beta*X(1,:,k)' );
%     lambda(1,:,k) = exp( mu_c + tempBeta*X(1,:,k)' );
end
%Y(1,:,:) = poissrnd(lambda(1,:,:));


% X = repmat(randn(T,size(X,2),1),[1 1 trials]);
% Run simulation
for k = 1:trials
    for i = 2:T
        % Determine change in X
        Noise(i,:,k) = mvnrnd(zeros(N,1),SigT);
%         dX(i,:,k) = permute( A*permute(X(i-1,:,k),[2 3 1])...
%             + Noise(i,:,k)', [3 1 2] )/tau;
        
        % Update X
%         X(i,:,k) = X(i-1,:,k) + dX(i,:,k);
        X(i,:,k) = (A*X(i-1,:,k)' + Noise(i,:,k)')/tau;
        
        % Generate spiking responses
        lambda(i,:,k) = exp( mu_c + Beta*X(i,:,k)' );
%         if i-BetaSize+1 < 2
%             tempBeta = reshape(Beta(end-i+1:end,:,:),[length(1:i)*N M])';
%             lambda(i,:,k) = exp( mu_c + ...
%                 tempBeta*reshape(X(1:i,:,k),numel(X(1:i,:,k)),1) + ...
%                 diag(Gamma(end-i+2:end,:)'*Y(1:i-1,:,k)) );
%         else
%             tempBeta = reshape(Beta,[BetaSize*N M])';
%             lambda(i,:,k) = exp( mu_c +...
%                 tempBeta*reshape(X(i-BetaSize+1:i,:,k),numel(X(i-BetaSize+1:i,:,k)),1) + ...
%                 diag(Gamma'*Y(i-GammaSize:i-1,:,k)) );
%         end

%         Y(i,:,k) = poissrnd(lambda(i,:,k));
        
    end
    for uniti = 1:M
        Y(:,uniti,k) = smooth(double(lambda(:,uniti,k) > rand(size(lambda(:,uniti,k)))),50);
    end
end
%Y = double(lambda > rand(size(lambda)));
% Y = poissrnd(lambda);

%% Fit model to the data
mu_ini = 0;
k_ini = [0 1 0 0 0 0 0 0 0 0];
h_ini = [-1 0 0 0 0 0 0 0 0 0];
optMethod = 'ML';
sigma = 0.5;
[mu, k, h, llike, exitflg, output] = LNP_fitter(X(:,:,:),Y(:,:,:),mu_ini,k_ini,h_ini,'optMethod',optMethod,'BasisSet',BasisSet);

%% Calculate some statistics
% Mean and variance of x(i)
Xbar = mean(X,3);
Xvar = var(X,[],3);

Ybar = mean(Y,3);%mean(Y,3);
Yvar = var(Y,[],3);%var(Y,[],3);

%% Perform PCA
Sigma = cov(Ybar);
[V, D] = eig(Sigma);
D = flipud(diag(D))/sum(diag(D));
V = fliplr(V);

YhatBar = Ybar*V;

for i = 1:trials
    Yhat(:,:,i) = Y(:,:,i)*V;
end
dYhat = diff(Yhat,1);

states = linspace(min(Yhat(:)),max(Yhat(:)),40);
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

% for indxi = 1:size(indx,1)
%     accept = BIN(:,1) == indx(indxi,1) & BIN(:,2) == indx(indxi,2);
%     n2(indxi,1) = sum(accept(:));
%     mdY(indxi,1) = mean(dYhat(accept,1));
%     mdY(indxi,2) = mean(dYhat(accept,2));
% end


%% Fit to the mean
Theta0 = [1 0 reshape(A,1,numel(A)) reshape(Beta,1,numel(Beta))]...
    + 0.000001*randn(size([1 0 reshape(A,1,numel(A)) reshape(Beta,1,numel(Beta))]));
constraintA = [zeros(1,N) reshape(diag(ones(1,N)),1,numel(diag(ones(1,N)))) zeros(1,N*M)];
constraintb = 0;        % Constrain trace(A) <= 0
lowerBound = [Theta0(1:N)-0.01 -ones(1,N^2)/10 -inf(1,N*M)];
upperBound = [Theta0(1:N)+0.01 ones(1,N^2)/10 inf(1,N*M)];
[x0, Ahat, Betahat, Xhat, Yhat] = fitLDSmean(Ybar,N,'Theta0',Theta0,...
    'constraintA',constraintA,'constraintb',constraintb,...
    'lowerBound',lowerBound,'upperBound',upperBound);

%% Plot the output

% Hidden States and measurements across trials
% figure('Name','Hidden states and measurements')
% trialinds = ceil(rand(2,1)*trials);
% for j = 1:N
%     subplot(N,1,j)
%     plot(squeeze(Y(:,j,trialinds)),'o','Color',[0.7 0.7 0.7])
%     hold on
%     plot(squeeze(X(:,j,trialinds)))
% end

% Plot a random subset of trials and the mean/var
figure('Name','Mean and var of X')
trialinds2 = ceil(rand(nexamps,1)*trials);
for j = 1:N
    subplot(N,1,j)
    h = myPatch(t,Xbar(:,j),sqrt(Xvar(:,j)),'patchProperties',patchProperties{j});
    hold on
    plot(t,squeeze(X(:,j,trialinds2)),'Color',colors(j,:))
    %plot(t,squeeze(Y(:,1,trialinds2(end))),'o','Color',colors(j,:))
    
    plot(t,Xbar(:,j),'Color',colors(j,:),'LineWidth',3)
end


figure('Name','Mean and var of Y')
trialinds2 = ceil(rand(nexamps,1)*trials);
for j = 1:3
    subplot(3,1,j)
    h = myPatch(t,Ybar(:,j),sqrt(Yvar(:,j)),'patchProperties',patchProperties{j});
    hold on
%     plot(t,squeeze(Y(:,j,trialinds2)),'Color',colors(j,:))
    
    plot(t,Ybar(:,j),'Color',colors(j,:),'LineWidth',3)
    axis tight
end

% PCA
figure('Name','PCA')
subplot(1,2,1)
contour(STATES{1},STATES{2},reshape(n2,size(STATES{1})))
colormap cool
hold on
allstates = [STATES{1}(:) STATES{2}(:)];
cutoff = 0.001;
quiver(allstates(n2/sum(n2) > cutoff,1),allstates(n2/sum(n2) > cutoff,2),...
    mdY(n2/sum(n2) > cutoff,1),mdY(n2/sum(n2) > cutoff,2),2,'Color',[0 0 0])
plot(YhatBar(:,1),YhatBar(:,2),'.-','LineWidth',1,'Color',[1 0 0])
plot(YhatBar(1,1),YhatBar(1,2),'o','MarkerSize',10,'Color',[1 0 0])
plot(YhatBar(end,1),YhatBar(end,2),'s','MarkerSize',10,'Color',[1 0 0])
axis([-3 3 -3 3])
axis square
mymakeaxis(gca)

subplot(1,2,2)
plot(cumsum(D),'ko')
axis square
mymakeaxis(gca)