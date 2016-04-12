%% SimulateHMM
%
% Simulates a linear-Gaussian hidden Markov model with N dimensions and M
% measurements.
%
%%

%% Parameters
% General
N = 1;              % Dimensionality of the system
M = 1;            % Measurement dimensionality
T = 1000;           % Total number of time-stamps
trials = 5000;       % Total number of trials
t = (1:T)';

% Transition model;
A = -1;%[0 -1; 1 0];        % Deterministic portion
SigT = 0.1;%[0.1 0; 0 0.1];      % Transition noise covariance
tau = 100;

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
mu_c = 0;
Beta = 0.01*flipud(BasisSet.k(2,:)');
Gamma = -0.005*flipud(BasisSet.h(1,:)');
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
spikes = nan(T,M,trials);        % Measurements;
X = nan(T,N,trials);        % Hidden state
dX = nan(T,N,trials);        % Change in hidden state

X(1,:,:) = repmat(rand(size(X(1,:,1))),[1 1 trials]);
tempBeta = reshape(Beta(end,:,:),[N, M])';
for k = 1:trials
    lambda(1,:,k) = exp( mu_c + tempBeta*X(1,:,k)' );
end
spikes(1,:,:) = poissrnd(lambda(1,:,:));


X = repmat(randn(T,size(X,2),1),[1 1 trials]);
% Run simulation
for k = 1:trials
    for i = 2:T
%         % Determine change in X
%         dX(i,:,k) = permute( A*permute(X(i-1,:,k),[2 3 1])...
%             + SigT*randn(N,1), [3 1 2] )/tau;
%         
%         % Update X
%         X(i,:,k) = X(i-1,:,k) + dX(i,:,k);
        
        % Generate spiking responses
        if i-BetaSize+1 < 2
            tempBeta = reshape(Beta(end-i+1:end,:,:),[length(1:i)*N M])';
            lambda(i,:,k) = exp( mu_c + ...
                tempBeta*reshape(X(1:i,:,k),numel(X(1:i,:,k)),1) + ...
                diag(Gamma(end-i+2:end,:)'*spikes(1:i-1,:,k)) );
        else
            tempBeta = reshape(Beta,[BetaSize*N M])';
            lambda(i,:,k) = exp( mu_c +...
                tempBeta*reshape(X(i-BetaSize+1:i,:,k),numel(X(i-BetaSize+1:i,:,k)),1) + ...
                diag(Gamma'*spikes(i-GammaSize:i-1,:,k)) );
        end
        spikes(i,:,k) = poissrnd(lambda(i,:,k));
        
    end
end

%% Fit model to the data
mu_ini = 0;
k_ini = [0 1 0 0 0 0 0 0 0 0];
h_ini = [-1 0 0 0 0 0 0 0 0 0];
optMethod = 'ML';
sigma = 0.5;
[mu, k, h, llike, exitflg, output] = LNP_fitter(X(:,:,:),spikes(:,:,:),mu_ini,k_ini,h_ini,'optMethod',optMethod,'BasisSet',BasisSet);

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