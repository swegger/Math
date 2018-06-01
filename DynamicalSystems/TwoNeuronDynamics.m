function TwoNeuronDynamics()
%% TwoNeuronDynamics
%
%   TwoNeuronDynamics()
%
%   Simulates the 2-neuron model from Wang et. al. 2017
%
%%

u0s = [0.6 0.6 0.65 0.65];
v0s = [0.65 0.65 0.6 0.6];
I0s = [0.6 0.8 0.6 0.8];

colors = [218 112 146;
           65 185 235;
          218 112 146;
           65 185 235]/255;

figure('Name','Two neuron','Position',[276 369 1231 495])
for i = 1:length(I0s)
    %% Parameters
    u0 = u0s(i);
    v0 = v0s(i);
    
    params.w = 6;
    params.epsi = 0.01;
    params.I0 = params.w * I0s(i);
    params.tau = 100;
    
    niter = 1001;
    
    vvals = linspace(0,1,100);
    
    %% Run simulation
    [u,v] = simulate_u_v(u0,v0,params,niter);
    
    
    %% Find null clines
    uvals = find_null_clines(vvals,params);
    
    %% Plot
    subplot(1,2,1)
    plot(u-v,'Color',colors(i,:))
    hold on
%     plot(v,'--','Color',colors(i,:))
    
    subplot(1,2,2)
    plot(u,v,'-','Color',colors(i,:))
    hold on
    plot(u(1:100:end),v(1:100:end),'.','Color',colors(i,:))
    plot(u(1),v(1),'ko','MarkerFaceColor','k','MarkerSize',7,'Color',colors(i,:))
    plot(u(end),v(end),'ks','MarkerFaceColor','k','MarkerSize',7,'Color',colors(i,:))
    
    plot(vvals,uvals(:,1),'Color',colors(i,:))
    plot(uvals(:,2),vvals,'--','Color',colors(i,:))
    axis([0 1.2 0 1.2])
    axis square

end

%% Clean up
subplot(1,2,1)
xlabel('t')
ylabel('u-v')
mymakeaxis(gca)

subplot(1,2,2)
plotUnity;
xlabel('u')
ylabel('v')
mymakeaxis(gca)

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
    
    uvals(:,1) = 1 ./ (1 + exp(-I0 + w*vvals - epsi) );
    uvals(:,2) = 1 ./ (1 + exp(-I0 + w*vvals ) );
    
    