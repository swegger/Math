%% SimulateHMM
%
% Simulates a linear-Gaussian hidden Markov model with N dimensions and M
% measurements.
%
%%

%% Parameters
% General
N = 2;              % Dimensionality of the system
M = 20;            % Measurement dimensionality
T = 1000;           % Total number of time-stamps
trials = 500;       % Total number of trials
t = (1:T)';

% Transition model;
A = [-0.000 -0.001; 0.5 -0.000];        % Deterministic portion
SigT = [0.001 0; 0 0.001];      % Transition noise covariance

% Measurement model
Beta = randn(M,N);
% Beta = zeros(M,N);
% Beta(1:floor(M/2),1) = 1;      % Meausrement gain
% Beta(floor(M/2)+1:end,2) = 1;      % Meausrement gain
SigM = 0.01.*diag(1:M);       % Measurement noise covariance

% Leak
Xl = 0;
gl = -0.1;

% Nonlinearity in X
f = @(x)( exp(x) );

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
Y = nan(T,M,trials);        % Measurements;
X = nan(T,N,trials);        % Hidden state

X(1,1:N,1:trials) = zeros(size(X(1,1:N,1:trials)));
X(1,1,1:trials) = ones(size(X(1,1,1:trials))) + 0.01*randn(size(X(1,1,1:trials)));
Y(1,1:M,1:trials) = permute( Beta*permute(X(1,1:N,1:trials),[2 3 1])...
    + SigM*randn(M,trials) , [3 1 2] );

% Run simulation
for i = 2:T
    % Update X
    X(i,1:N,1:trials) = X(i-1,1:N,1:trials)...
        + gl*(X(i-1,1:N,1:trials) - Xl)...
        + permute( A*permute(f(X(i-1,1:N,1:trials)),[2 3 1])...
        + SigT*randn(N,trials), [3 1 2] );
    
    % Generate measurements
    Y(i,1:M,1:trials) = permute( Beta*permute(X(i,1:N,1:trials),[2 3 1])...
        + SigM*randn(M,trials) , [3 1 2] );
    
end

%% Calculate some statistics
% Mean and variance of x(i)
Xbar = mean(X,3);
Xvar = var(X,[],3);

Ybar = mean(Y,3);
Yvar = var(Y,[],3);

%% Perform PCA
Sigma = cov(Ybar');
[V, D] = eig(Sigma);
D = flipud(diag(D))/sum(diag(D));
V = fliplr(V);

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
figure('Name','Hidden states and measurements')
trialinds = ceil(rand(2,1)*trials);
for j = 1:N
    subplot(N,1,j)
    plot(squeeze(Y(:,j,trialinds)),'o','Color',[0.7 0.7 0.7])
    hold on
    plot(squeeze(X(:,j,trialinds)))
end

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
for j = 1:M
    subplot(M,1,j)
    h = myPatch(t,Ybar(:,j),sqrt(Yvar(:,j)),'patchProperties',patchProperties{j});
    hold on
    plot(t,squeeze(Y(:,j,trialinds2)),'Color',colors(j,:))
    %plot(t,squeeze(Y(:,1,trialinds2(end))),'o','Color',colors(j,:))
    
    plot(t,Ybar(:,j),'Color',colors(j,:),'LineWidth',3)
end