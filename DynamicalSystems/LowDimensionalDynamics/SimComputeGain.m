function SimComputeGain(varargin)
%%
%
%
%%

%% Defaults
K_default = [0 0 0; 1 1 0; 0.7 0.3 0];    % Defaults to a low pass of input and an integrator of the low-pass
f_default = @(v)(v); %@(v)(1./(1+exp(-v)))

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'N',100)
addParameter(Parser,'K',K_default)
addParameter(Parser,'Gamma',NaN)
addParameter(Parser,'zeta',NaN)
addParameter(Parser,'t',-100:1600)
addParameter(Parser,'trialN',50)
addParameter(Parser,'Sigma',NaN)
addParameter(Parser,'f',f_default)
addParameter(Parser,'tau',300)
addParameter(Parser,'constant',0)

parse(Parser,varargin{:})

N = Parser.Results.N;
K = Parser.Results.K;
Gamma = Parser.Results.Gamma;
zeta = Parser.Results.zeta;
t = Parser.Results.t;
trialN = Parser.Results.trialN;
Sigma = Parser.Results.Sigma;
f = Parser.Results.f;
tau = Parser.Results.tau;
constant = Parser.Results.constant;

M = size(K,1);

% Generate compression matrices
if any(isnan(Gamma(:)))
    Gamma = randn(M,N)/sqrt(N);
end

% Generate input
if any(isnan(zeta(:)))
    zeta = generateZeta(M,t);
end

% Set covariance matrix of inputs
if any(isnan(Sigma(:)))
    sigma = 0.1;
    Sigma = diag(sigma*ones(N,1));
end

deltat = t(2)-t(1);

%% Simulate voltages and spike counts
W = Gamma'*(K-eye(M))*Gamma;
u = Gamma'*zeta;

v = nan(N,length(t),trialN);
dv = nan(N,length(t),trialN);
nu = nan(M,length(t),trialN);

v(:,1,:) = randn([N,1,trialN])*sigma;
dv(:,1,:) = zeros([N,1,trialN]);


v2 = nan(N,length(t),trialN);
dv2 = nan(N,length(t),trialN);
nu2 = nan(M,length(t));
dnu2 = nan(M,length(t));

v2(:,1,:) = randn([N,1,trialN])/10;
dv2(:,1,:) = zeros([N,1,trialN]);
nu2(:,1,:) = zeros([M,1]);
dnu2(:,1,:) = zeros([M,1]);

for ti = 2:length(t)
    dnu2(:,ti) = ...
        ((K-eye(M))*nu2(:,ti-1) + zeta(:,ti))/tau;
    nu2(:,ti) = nu2(:,ti-1) + dnu2(:,ti-1)*deltat;
end

for triali = 1:trialN
    eta = permute(mvnrnd(zeros(length(t),N),Sigma),[2,1]);
    for ti = 2:length(t)
        dv(:,ti,triali) = (...
            W*f(v(:,ti-1,triali)) + u(:,ti) + eta(:,ti) + constant)/tau;
        v(:,ti,triali) = v(:,ti-1,triali) + dv(:,ti,triali)*deltat;
        
        dv2(:,ti,triali) = (-v2(:,ti-1,triali) + ...
            Gamma'*nu2(:,ti) + u(:,ti)/tau + eta(:,ti)/tau + constant/tau);
        v2(:,ti,triali) = v2(:,ti-1,triali) + dv2(:,ti,triali)*deltat;
    end
    nu(:,:,triali) = Gamma*f(v(:,:,triali));
end

r = f(v);

%% Plotting
figure
for mi = 1:M
    subplot(1,M,mi)
    plot(t,permute(nu(mi,:,randsample(trialN,10)),[2,3,1]),'Color',[0.6 0.6 0.6])
    hold on
    plot(t,permute(mean(nu(mi,:,:),3),[2,1]),'k','LineWidth',2)
    plot(t,permute(nu2(mi,:),[2,1]),'r','LineWidth',1)
end

%% functions

%% generateZeta
function zeta = generateZeta(M,t)
    zeta = zeros(M,length(t));
    zeta(1,t >= 0 & t < 600) = 1;
    
%% localSigmoid
function r = localSigmoid(v)
    r = 1./(1+exp(v));