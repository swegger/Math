%% exponentialModels
%
%   Simulates firing rate models of the form:
%       r(t) = f( (v(t)-vT)/vS )
%       dv/dt = -(v(t)-El) + G*r(t) + eta(t);
%
%%  

%% Variables

% General
N = 2;              % Dimensionality of the system
T = 500;           % Total number of time-stamps
trials = 1;       % Total number of trials
t = (1:T)';

% Coupling terms
G = [0 1; -1 0];
%g = 4;
%sigG = g/sqrt(N);
%G = sigG*randn(N,N);
%G = G - diag(diag(G));

% Nonlinear terms
vT = 1;         % (soft) Threshold
vT2 = 2;        % (soft) Threshold 2
vS = 1;         % Sensitivity

% Time constant
tau = 10;

% Leak terms
gl = 0;
El = 0;

% Constant input terms
I0 = 0;
I02 = 0;

% Noise Terms
SigT = diag(0.1*ones(1,N)) * 0;

% Firing rate function
%f= @(x,T,S)(x);
%f = @(x,T,S)( tanh( (x-T)/S ) );
f = @(x,T,S)( exp( (x - T)/S ) );
spkf = @(lambda)( lambda );
%spkf = @(lambda)( poissrnd(lambda) );
spkmax = Inf;

% Gain function
g = @(u)(u);
u = 1;
u2 = 1;

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
dv = nan(T,N,trials);
v = nan(T,N,trials);       
r = nan(T,N,trials);      
r2 = nan(T,N,trials); 
spk = nan(T,N,trials);
spk2 = nan(T,N,trials);

dv(1,:,:) = zeros(size(dv(1,:,:)));
dv2 = dv;
% v(1,:,:) = repmat(2*rand(size(v(1,:,1))) - 1,[1 1 trials]);
v(1,:,:) = repmat(zeros(size(v(1,:,1))),[1 1 trials]);
v2 = v;
r(1,:,:) = f(v(1,:,:),vT,vS);
r2(1,:,:) = f(v2(1,:,:),vT2,vS);
spk(1,:,:) = spkf(r(1,:,:));
spk2(1,:,:) = spkf(r2(1,:,:));
spk(spk>spkmax) = spkmax;
spk2(spk>spkmax) = spkmax;

% Simulate the model
for i = 2:T
    % Update voltages
    dv(i,1:N,1:trials) = ( -gl*(v(i-1,1:N,1:trials) - El)...
        + permute( G*permute(spk(i-1,1:N,1:trials),[2 3 1])...
        + SigT*randn(N,trials), [3 1 2] )...
        + I0 )/tau;
    v(i,1:N,1:trials) = v(i-1,1:N,1:trials)...
        + dv(i,1:N,1:trials);
    
    dv2(i,1:N,1:trials) = ( -gl*(v2(i-1,1:N,1:trials) - El)...
        + permute( G*permute(spk2(i-1,1:N,1:trials),[2 3 1])...
        + SigT*randn(N,trials), [3 1 2] )...
        + I02 )/tau;
    v2(i,1:N,1:trials) = v2(i-1,1:N,1:trials)...
        + dv2(i,1:N,1:trials);
    
    % Generate firing rates
    r(i,:,:) = g(u)*f(v(i,:,:),vT,vS);
    r2(i,:,:) = g(u2)*f(v2(i,:,:),vT2,vS);
    spk(i,:,:) = spkf(r(i,:,:));
    spk(spk>spkmax) = spkmax;
    spk2(i,:,:) = spkf(r2(i,:,:));
    spk2(spk>spkmax) = spkmax;
end

%% Calculate statistics
% Mean and variance of x(i)
vbar = nanmean(v,3);
vbar2 = nanmean(v2,3);
vvar = nanvar(v,[],3);
vvar2 = nanvar(v2,[],3);

rbar = nanmean(r,3);
rbar2 = nanmean(r2,3);
rvar = nanvar(r,[],3);
rvar2 = nanvar(r2,[],3);

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
unitinds = [1 2]; %ceil(rand(nexamps,1)*N);
trialinds2 = ceil(rand(nexamps,1)*trials);
for j = 1:length(unitinds)
    subplot(length(unitinds),1,j)
    h = myPatch(t,vbar(:,unitinds(j)),sqrt(vvar(:,unitinds(j))),'patchProperties',patchProperties{j});
    hold on
    plot(t,squeeze(v(:,unitinds(j),trialinds2)),'Color',colors(j,:))
    
    plot(t,vbar(:,unitinds(j)),'Color',colors(j,:),'LineWidth',3)
    plot(t,vbar2(:,unitinds(j)),'--','Color',colors(j,:),'LineWidth',3)
    ax = axis;
    %axis([ax(1:2) -5 5])
end


figure('Name','Mean and var of r')
for j = 1:length(unitinds)
    subplot(length(unitinds),1,j)
    %h = myPatch(t,rbar(:,unitinds(j)),sqrt(rvar(:,unitinds(j))),'patchProperties',patchProperties{j});
    hold on
    plot(t,squeeze(r(:,unitinds(j),trialinds2)),'Color',colors(j,:))
    
    plot(t,rbar(:,unitinds(j)),'Color',colors(j,:),'LineWidth',3)
    plot(t,f(vbar2(:,unitinds(j)),vT,vS),'--','Color',colors(j,:),'LineWidth',3)
    xlabel('Time')
    ylabel('r')
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

% Plot the f-I curve
figure('Name','f-I')
plot(v(:),r(:),'.')
hold on
plot(v2(:),r2(:),'.')
