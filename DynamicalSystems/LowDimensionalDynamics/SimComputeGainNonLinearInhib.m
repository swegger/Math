function SimComputeGainNonLinearInhib(varargin)
%%
%
%
%%

%% Defaults
K_default = [0 0 0 0; 1 1 0 0; 0.7 0.3 0 -0.1; 0 0 1 0];    % Defaults to a low pass of input and an integrator of the low-pass
f_default = @(v)(v); %@(v)(1./(1+exp(-v)))
f_coherence_default = @(c)(1-c/3);

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'N',100)
addParameter(Parser,'K',K_default)
addParameter(Parser,'Gamma',NaN)
addParameter(Parser,'zeta',NaN)
addParameter(Parser,'t',-100:1600)
addParameter(Parser,'trialN',50)
addParameter(Parser,'cohs',[20 60 100]/100)
addParameter(Parser,'s',10)
addParameter(Parser,'Sigma',NaN)
addParameter(Parser,'f',f_default)
addParameter(Parser,'f_coherence',f_coherence_default)
addParameter(Parser,'tau',300)
addParameter(Parser,'betas',[0.4 -0.1])
addParameter(Parser,'constant',0)
addParameter(Parser,'baseline',log(10))
addParameter(Parser,'sparsity',NaN)
addParameter(Parser,'inhibInd',4)

parse(Parser,varargin{:})

N = Parser.Results.N;
K = Parser.Results.K;
Gamma = Parser.Results.Gamma;
zeta = Parser.Results.zeta;
t = Parser.Results.t;
trialN = Parser.Results.trialN;
cohs = Parser.Results.cohs;
s = Parser.Results.s;
Sigma = Parser.Results.Sigma;
f = Parser.Results.f;
f_coherence = Parser.Results.f_coherence;
tau = Parser.Results.tau;
betas = Parser.Results.betas;
constant = Parser.Results.constant;
baseline = Parser.Results.baseline;
sparsity = Parser.Results.sparsity;
inhibInd = Parser.Results.inhibInd;

M = size(K,1);

coherence = randsample(cohs,trialN,true);
speeds = randsample(s,trialN,true);

% Generate compression matrices
if any(isnan(Gamma(:)))
    Gamma = randn(M,N)/sqrt(N);
end
Gamma(abs(Gamma)<sparsity) = 0;

% Generate input
if any(isnan(zeta(:)))
    zeta = generateZeta(M,t,coherence,betas);
end

% Set covariance matrix of inputs
if any(isnan(Sigma(:)))
    sigma = 0.1;
    Sigma = diag(sigma*ones(N,1));
end

deltat = t(2)-t(1);

%% Simulate voltages and spike counts
W = Gamma'*(K-eye(M))*Gamma;
for triali = 1:trialN
    u(:,:,triali) = Gamma'*zeta(:,:,triali);
end
v = nan(N,length(t),trialN);
nu = nan(M,length(t),trialN);

v(:,1,:) = randn([N,1,trialN])*sigma;


v2 = nan(N,length(t),trialN);
nu2 = nan(M,length(t),trialN);

v2(:,1,:) = randn([N,1,trialN])/10;
nu2(:,1,:) = zeros([M,1,trialN]);

for triali = 1:trialN
    eta = permute(mvnrnd(zeros(length(t),N),Sigma),[2,1]);
    for ti = 2:length(t)
        coherence_mod = zeros(M,1);
        coherence_mod(inhibInd) = f_coherence(coherence(triali));
%         dv = (...
%             W*f(v(:,ti-1,triali)+Gamma'*coherence_mod) + u(:,ti,triali) + eta(:,ti) + constant)/tau;
        x = Gamma*v(:,ti-1,triali);
        dx = (K-eye(M))*x;
        dx(inhibInd) = -x(inhibInd) + K(inhibInd,:)*x*f_coherence(coherence(triali));
        dv = (Gamma'*dx + u(:,ti,triali) + eta(:,ti) + constant)/tau;
        v(:,ti,triali) = v(:,ti-1,triali) + dv*deltat;
        
        dnu2 = ...
            ((K-eye(M))*nu2(:,ti-1,triali) + zeta(:,ti,triali))/tau;
        nu2(:,ti,triali) = nu2(:,ti-1,triali) + dnu2*deltat;
        
        dv2 = (-v2(:,ti-1,triali) + ...
            Gamma'*nu2(:,ti,triali) + u(:,ti,triali)/tau + eta(:,ti)/tau + constant/tau);
        v2(:,ti,triali) = v2(:,ti-1,triali) + dv2*deltat;
    end
    nu(:,:,triali) = Gamma*f(v(:,:,triali));
end

counts = poissrnd(exp(v + baseline));
countsCI = poissrnd(exp(v2 + baseline));

%% Find mean and covariance of counts for each neuron
cohs = unique(coherence);
sz = size(counts);
res = nan(sz(2),sz(1),sz(3));
resCI = nan(sz(2),sz(1),sz(3));
for ni = 1:N
    ind = 1;
    for ci = 1:length(cohs)
        m(:,ni,ci) = mean(counts(ni,:,coherence == cohs(ci)),3);
        va(:,ni,ci) = var(counts(ni,:,coherence == cohs(ci)),[],3);
        mCI(:,ni,ci) = mean(countsCI(ni,:,coherence == cohs(ci)),3);
        vaCI(:,ni,ci) = var(countsCI(ni,:,coherence == cohs(ci)),[],3);
        
        resTemp = permute(counts(ni,:,coherence == cohs(ci)),[2,1,3]) - ...
            repmat(m(:,ni,ci),[1,1,sum(coherence == cohs(ci))]);
        res(:,ni,ind:ind+sum(coherence == cohs(ci))-1) = resTemp;
        
        
        resTemp = permute(countsCI(ni,:,coherence == cohs(ci)),[2,1,3]) - ...
            repmat(m(:,ni,ci),[1,1,sum(coherence == cohs(ci))]);
        resCI(:,ni,ind:ind+sum(coherence == cohs(ci))-1) = resTemp;
        
        ind = ind+sum(coherence == cohs(ci));
        
        C(:,:,ni,ci) = cov(permute(counts(ni,:,coherence == cohs(ci)),[3,2,1]));
    end
end

Me = permute(mean(counts,3),[2,1]);
V = var(res,[],3);
FF = V./Me;
phi = min(FF,[],1);
varCE = V - repmat(phi,[size(Me,1),1]).*Me;

MeCI = permute(mean(countsCI,3),[2,1]);
VCI = var(resCI,[],3);
FFCI = VCI./MeCI;
phiCI = min(FFCI,[],1);
varCE_CI = VCI - repmat(phiCI,[size(MeCI,1),1]).*MeCI;

for ni = 1:N
    Cres(:,:,ni) = cov(permute(res(:,ni,:),[3,1,2]));
    CresCI(:,:,ni) = cov(permute(resCI(:,ni,:),[3,1,2]));
end

%% PCA
CR = cov(reshape( permute(m - mean(m,[1,3]),[1,3,2]),[size(m,1)*size(m,3),size(m,2)]));
[vecs, vals] = eig(CR);
vecs = fliplr(vecs);
vals = flipud(diag(vals));
reconstructionN = 5;

for ci = 1:size(m,3)
    recon(:,:,ci) = permute(vecs(:,1:reconstructionN)'*m(:,:,ci)',[2,3,4,1]);
end

%% Fitting
ind = 19;
x = reshape(nu2,[size(nu2,1),size(nu2,2)*size(nu2,3)]);
y = reshape(counts(ind,:,:),[1,size(counts,2)*size(counts,3)]);
[theta, logPosterior, exitflg, output, thetas] = FitMultivariateLinearExpPoissonObsMAP(x,y,[0 0 0 0]);

%% Plotting
figure('Name','Low-dimensional computations','Position',[336 503 1754 420])
colors = colormap('lines');
for mi = 1:M
    subplot(1,M,mi)
    plot(t,permute(nu(mi,:,randsample(trialN,10)),[2,3,1]),'Color',[0.6 0.6 0.6])
    hold on
    
    for ci = 1:length(cohs)
        plot(t,permute(mean(nu(mi,:,coherence == cohs(ci)),3),[2,1]),...
            'Color',colors(ci,:),'LineWidth',2)
    end
    xlabel('Time form input onset')
    ylabel('Low-dimensional output')
end

figure('Name','Computational covariance','Position',[1262 367 1180 962])
subplot(2,2,1)
imagesc(t,t,mean(Cres,3)-diag(diag(mean(Cres,3))))%+diag(mean(varCE,2)))
cax = caxis;
hold on
axis square
plotVertical([0 600]);
plotHorizontal([0 600]);
xlabel('Time from input onset')
ylabel('Time from input onset')
title('CovCE (network dynamics)')

subplot(2,2,2)
imagesc(t,t,mean(CresCI,3)-diag(diag(mean(CresCI,3))))%+diag(mean(varCE_CI,2)))
cax2 = caxis;
hold on
axis square
plotVertical([0 600]);
plotHorizontal([0 600]);
xlabel('Time from input onset')
ylabel('Time from input onset')
title('CovCE (common input)')
caxis([min(cax(1),cax2(1)) max(cax(2),cax2(2))])

subplot(2,2,1)
caxis([min(cax(1),cax2(1)) max(cax(2),cax2(2))])

subplot(2,2,[3 4])
plot(t,mean(varCE,2))
hold on
plot(t,mean(varCE_CI,2))
axis tight
plotVertical([0 600]);
xlabel('Time from input onset')
ylabel('varCE')
legend({'Network','Common Input'},'Location','Northwest')


%% functions

%% generateZeta
function zeta = generateZeta(M,t,coherence,betas)
    zeta = zeros(M,length(t),length(coherence));
%     utemp = repmat(permute(coherence,[1,3,2]),[1,sum(t>=0 & t<600),1]);
%     zeta(1,t >= 0 & t < 600,:) = utemp;
    
    zetaBase = zeta;
    zetaBase(1,t >= 0 & t < 200,:) = ones(size(zetaBase(1,t >= 0 & t < 200,:)));
    zetaBase(1,t >= 200,:) = 0.1*ones(size(zetaBase(1,t >= 200,:)));
    for ci = 1:length(coherence)
        zeta(1,t >= 0 & t < 200,ci) = zetaBase(1,t >= 0 & t < 200,ci) + betas(1)*(coherence(ci)-mean(coherence));
        zeta(1,t >= 200,ci) = zetaBase(1,t >= 200,ci) + betas(2)*(coherence(ci)-mean(coherence));
    end
    
%% localSigmoid
function r = localSigmoid(v)
    r = 1./(1+exp(v));
    
%% generateK
function K = generateK(Kbase,cohAdjust,inds)
    K = Kbase;
    for ii = size(inds,1)
        K(inds(ii,1),inds(ii,2)) = cohAdjust(ii);
    end
    