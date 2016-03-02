%% exponentialModels
%
%   Simulates firing rate models of the form:
%       r(t) = f( (v(t)-vT)/vS )
%       dv/dt = -(v(t)-El) + G*r(t) + eta(t);
%
%%  

%% Variables

% General
Ne = 500;              % Dimensionality of the system
Ni = 500;              % Dimensionality of the system
T = 1000;           % Total number of time-stamps
trials = 1;       % Total number of trials
t = (1:T)';

% Coupling terms
g = 1.5;
sigGe = g/sqrt(Ne);
sigGi = g/sqrt(Ni);
Gee = sigGe*abs(randn(Ne,Ne));
Gei = -sigGi*abs(randn(Ne,Ni));
Gie = sigGe*abs(randn(Ni,Ne));
Gii = -sigGi*abs(randn(Ni,Ni));
% Gee = 0.0;
% Gei = -0.1;
% Gie = 0.05;
% Gii = -0;

% Time constants
tau = 100;

% Nonlinear terms
vT = 0;         % (soft) Threshold
vT2 = 0;        % (soft) Threshold 2
vS = 1;         % Sensitivity
vS2 = 1;         % Sensitivity

% Leak terms
gl = 1;
El = 0;

% Constant input terms
I0 = 0;
I02 = 0;

% Noise Terms
SigTe = diag(0.1*ones(1,Ne)) * 0;
SigTi = diag(0.1*ones(1,Ni)) * 0;

% Firing rate function
%f= @(x,T,S)(x);
f = @(x,T,S)( tanh( (x-T)/S ) );
%f = @(x,T,S)( exp( (x - T)/S ) );
spkf = @(lambda)( lambda );
spkmax = Inf;
%spkf = @(lambda)( poissrnd(lambda) );

% Plot options
nexamps = 5;    % Number of example trials to show in mean and variance plot
Neigs = 3;
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
for i = 1:size(colors,1)
    patchProperties{i}.FaceColor = colors(i,:) + 0.5*(ones(1,3)-colors(i,:));
end

%% Run the simulation

% Initialize the matrices
ve = nan(T,Ne,trials);       
re = nan(T,Ne,trials);      
re2 = nan(T,Ne,trials); 
spke = nan(T,Ne,trials);
spke2 = nan(T,Ne,trials);

vi = nan(T,Ni,trials);       
ri = nan(T,Ni,trials);      
ri2 = nan(T,Ni,trials); 
spki = nan(T,Ni,trials);
spki2 = nan(T,Ni,trials);

ve(1,:,:) = repmat(2*rand(size(ve(1,:,1))) - 1,[1 1 trials]);
ve2 = ve;
re(1,:,:) = f(ve(1,:,:),vT,vS);
re2(1,:,:) = f(ve(1,:,:),vT2,vS2);
spke(1,:,:) = spkf(re(1,:,:));
spke2(1,:,:) = spkf(re2(1,:,:));
spke(spke>spkmax) = spkmax;
spke2(spke>spkmax) = spkmax;

vi(1,:,:) = repmat(2*rand(size(vi(1,:,1))) - 1,[1 1 trials]);
vi2 = vi;
ri(1,:,:) = f(vi(1,:,:),vT,vS);
ri2(1,:,:) = f(vi(1,:,:),vT2,vS2);
spki(1,:,:) = spkf(ri(1,:,:));
spki2(1,:,:) = spkf(ri2(1,:,:));
spki(spki >spkmax) = spkmax;
spki2(spki>spkmax) = spkmax;

% Simulate the model
for i = 2:T
    % Update voltages
    ve(i,1:Ne,1:trials) = ve(i-1,1:Ne,1:trials)...
        + ( -gl*(ve(i-1,1:Ne,1:trials) - El)...
        + permute( Gee*permute(spke(i-1,1:Ne,1:trials),[2 3 1])...
        + Gei*permute(spki(i-1,1:Ni,1:trials),[2 3 1])...
        + SigTe*randn(Ne,trials), [3 1 2] )...
        + I0 )/tau;
    
    ve2(i,1:Ne,1:trials) = ve2(i-1,1:Ne,1:trials)...
        + ( -gl*(ve2(i-1,1:Ne,1:trials) - El)...
        + permute( Gee*permute(spke2(i-1,1:Ne,1:trials),[2 3 1])...
        + Gei*permute(spki2(i-1,1:Ni,1:trials),[2 3 1])...
        + SigTe*randn(Ne,trials), [3 1 2] )...
        + I02 )/tau;
    
    vi(i,1:Ni,1:trials) = vi(i-1,1:Ni,1:trials)...
        + ( -gl*(vi(i-1,1:Ni,1:trials) - El)...
        + permute( Gii*permute(spki(i-1,1:Ni,1:trials),[2 3 1])...
        + Gie*permute(spke(i-1,1:Ne,1:trials),[2 3 1])...
        + SigTi*randn(Ni,trials), [3 1 2] )...
        + I0 )/tau;
    
    vi2(i,1:Ni,1:trials) = vi2(i-1,1:Ni,1:trials)...
        + ( -gl*(vi2(i-1,1:Ni,1:trials) - El)...
        + permute( Gii*permute(spki2(i-1,1:Ni,1:trials),[2 3 1])...
        + Gie*permute(spke2(i-1,1:Ne,1:trials),[2 3 1])...
        + SigTi*randn(Ni,trials), [3 1 2] )...
        + I02  )/tau;
    
    % Generate firing rates
    re(i,:,:) = f(ve(i,:,:),vT,vS);
    re2(i,:,:) = f(ve2(i,:,:),vT2,vS2);
    spke(i,:,:) = spkf(re(i,:,:));
    spke(spke>spkmax) = spkmax;
    spke2(i,:,:) = spkf(re2(i,:,:));
    spke2(spke>spkmax) = spkmax;
    
    ri(i,:,:) = f(vi(i,:,:),vT,vS);
    ri2(i,:,:) = f(vi2(i,:,:),vT2,vS2);
    spki(i,:,:) = spkf(ri(i,:,:));
    spki(spki>spkmax) = spkmax;
    spki2(i,:,:) = spkf(ri2(i,:,:));
    spki2(spki>spkmax) = spkmax;
end

%% Calculate statistics
% Mean and variance of x(i)
vbar = nanmean(ve,3);
vbar2 = nanmean(ve2,3);
vvar = nanvar(ve,[],3);
vvar2 = nanvar(ve2,[],3);

rbar = nanmean(re,3);
rbar2 = nanmean(re2,3);
rvar = nanvar(re,[],3);
rvar2 = nanvar(re2,[],3);

%% Perform PCA
Sigma = nancov(rbar');
[V, D] = eig(Sigma);
D = flipud(diag(D))/sum(diag(D));
V = fliplr(V);

Sigma = nancov(rbar2');
[V2, D2] = eig(Sigma);
D2 = flipud(diag(D2))/sum(diag(D2));
V2 = fliplr(V2);

%% Plot the output

% Plot a random subset of trials and the mean/var
figure('Name','Mean and var of v')
unitinds = ceil(rand(nexamps,1)*Ne);
trialinds2 = ceil(rand(nexamps,1)*trials);
for j = 1:length(unitinds)
    subplot(length(unitinds),1,j)
    h = myPatch(t,vbar(:,unitinds(j)),sqrt(vvar(:,unitinds(j))),'patchProperties',patchProperties{j});
    hold on
    plot(t,squeeze(ve(:,unitinds(j),trialinds2)),'Color',colors(j,:))
    
    plot(t,vbar(:,unitinds(j)),'Color',colors(j,:),'LineWidth',3)
    plot(t,vbar2(:,unitinds(j)),'--','Color',colors(j,:),'LineWidth',3)
    ax = axis;
    %axis([ax(1:2) -5 5])
end


figure('Name','Mean and var of r')
for j = 1:length(unitinds)
    subplot(length(unitinds),1,j)
    h = myPatch(t,rbar(:,unitinds(j)),sqrt(rvar(:,unitinds(j))),'patchProperties',patchProperties{j});
    hold on
    plot(t,squeeze(re(:,unitinds(j),trialinds2)),'Color',colors(j,:))
    
    plot(t,rbar(:,unitinds(j)),'Color',colors(j,:),'LineWidth',3)
    plot(t,rbar2(:,unitinds(j)),'--','Color',colors(j,:),'LineWidth',3)
end

% Plot the top 3 eigenvectors for the different thresholds
figure('Name','Eigenvectors')
for j = 1:Neigs
    subplot(Neigs,1,j)
    plot(t,V(:,j),'k','LineWidth',3)
    hold on
    plot(t,V2(:,j),'k--','LineWidth',3)
    ax(j,:) = axis;
end
ax = [min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))];
for j = 1:Neigs
    subplot(Neigs,1,j)
    axis(ax);
    text(T/2,ax(3)/2,['- % Var = ' num2str(D(j))])
    text(T/2,ax(3)/1.5,['-- % Var = ' num2str(D2(j))])
end