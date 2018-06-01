function TwoNeuronDynamics()
%% TwoNeuronDynamics
%
%   TwoNeuronDynamics()
%
%   Simulates the 2-neuron model from Wang et. al. 2017
%
%%

threshold = 0.4;
u0s = [0.6 0.6 0.65 0.65];
v0s = [0.65 0.65 0.6 0.6];
I0s = [0.6 0.8 0.6 0.8];

colors = [218 112 146;
           65 185 235;
          218 112 146;
           65 185 235]/255;

h1 = figure('Name','Two neuron','Position',[204 160 528 918]);
h2 = figure('Name','Two neuron','Position',[204 160 999 918]);
for i = 1:length(I0s)
    %% Parameters
    u0 = u0s(i);
    v0 = v0s(i);
    
    params.w = 6;
    params.epsi = 0.01;
    params.I0 = params.w * I0s(i);
    params.tau = 100;
    
    niter = 1001;
    
    vvals = linspace(0,1.2,100);
    
    %% Run simulation
    [u,v] = simulate_u_v(u0,v0,params,niter);
    
    
    %% Find null clines
    uvals = find_null_clines(vvals,params);
    
    %% Plot
    %subplot(3,3,[1 4 7])
    figure(h1)
%     if i >= 3
        temp = u-v;
        plot(temp(abs(temp) < threshold),'Color',colors(i,:))
        hold on
%     end
%     plot(v,'--','Color',colors(i,:))
    
    figure(h2)
    if i < 3
        subplot(1,2,1)
    else
        subplot(1,2,2)
    end
%    subplot(ceil(sqrt(length(I0s))),ceil(sqrt(length(I0s))),i)
%     plot(u,v,'-','Color',colors(i,:))
    hold on
    plot(u(1:50:end),v(1:50:end),'.','Color',colors(i,:))
    plot(u(1),v(1),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',7,'Color',colors(i,:))
    plot(u(end),v(end),'ko','MarkerFaceColor',colors(i,:),'MarkerSize',7,'Color','k')
    
    plot(vvals,uvals(:,1),'Color',colors(i,:))
    plot(uvals(:,2),vvals,'--','Color',colors(i,:))
    
    axis([0 1.2 0 1.2])
    axis square
    plotUnity;
    xlabel('u')
    ylabel('v')
    mymakeaxis(gca)


end

figure(h1)
xlabel('Time')
ylabel('x_1-x_2')
mymakeaxis(gca,'yticks',[-1 0 1])

%% Functions

%% Activation function
function out = thresh_exp(x)
    out = 1 ./ (1 + exp(-x));
    
%% find_u_dot
function du = find_u_dot(u,v,params)
    w = params.w;
    epsi = params.epsi;
    I0 = params.I0;
    tau = params.tau;
    
    du = ( -u + thresh_exp( I0 - w*v ) )/tau;
    
%% find_v_dot
function dv = find_v_dot(u,v,params)
    w = params.w;
    epsi = params.epsi;
    I0 = params.I0;
    tau = params.tau;
    
    dv = ( -v + thresh_exp( I0 - w*u ) )/tau;

%% simulate_u_v
function [u, v] = simulate_u_v(u0,v0,params,niter)
    
    u = nan(1,niter);
    v = nan(1,niter);
    
    u(1) = u0;
    v(1) = v0;
    for i = 2:niter
        u(i) = u(i-1) + find_u_dot(u(i-1),v(i-1),params);
        v(i) = v(i-1) + find_v_dot(u(i-1),v(i-1),params);
    end
    
%% find_null_clines
function uvals = find_null_clines(vvals,params)
    w = params.w;
    I0 = params.I0;
    epsi = params.epsi;
    
    uvals(:,1) = 1 ./ (1 + exp(-I0 + w*vvals ) );
    uvals(:,2) = 1 ./ (1 + exp(-I0 + w*vvals ) );
    
    