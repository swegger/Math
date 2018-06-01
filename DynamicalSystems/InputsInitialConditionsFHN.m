function InputsInitialConditionsFHN()
%% InputsInitialConditionsFHN
%
%   Simulates Fitzhugh-Nagomo model (a al Ringal 2006) for 4 different
%   initial conditions and 
%%

% Constant variables
dt = 0.01;              % Time step of simulation
epsilon = 0.01;         % Relative time scales of each dimension
b = [1 0.1];            % 'weight matrix'
uvals = linspace(-0.5,1.2,20);   % For vector fields and nullclines
wvals = linspace(-0.15,0.5,20);

t = 0:dt:400;

% IC(1), I(1)
vars(1).I = 0;                  % External input
vars(1).u0 = -0.1;                % Inital condition of u
vars(1).w0 = 0;                % Initial condition of w

% IC(2), I(1);
vars(2).I = 0;                  % External input
vars(2).u0 = -0.1;                % Inital condition of u
vars(2).w0 = -0.1;                % Initial condition of w

% IC(1), I(2)
vars(3).I = 0.075;                  % External input
vars(3).u0 = -0.1;                % Inital condition of u
vars(3).w0 = 0;                % Initial condition of w

% IC(1), I(2)
vars(4).I = 0.2;                  % External input
vars(4).u0 = -0.1;                % Inital condition of u
vars(4).w0 = 0;                % Initial condition of w

% % IC(2), I(2);
% vars(4).I = 0.1;                  % External input
% vars(4).u0 = -0.1;                % Inital condition of u
% vars(4).w0 = -0.1;                % Initial condition of w

% % IC(1), I(1)
% vars(1).I = 0;                  % External input
% vars(1).u0 = -0.1;                % Inital condition of u
% vars(1).w0 = -0.01;                % Initial condition of w
% 
% % IC(2), I(1);
% vars(2).I = 0;                  % External input
% vars(2).u0 = -0.1;                % Inital condition of u
% vars(2).w0 = -0.1;                % Initial condition of w
% 
% % IC(1), I(2)
% vars(3).I = 0.05;                  % External input
% vars(3).u0 = -0.1;                % Inital condition of u
% vars(3).w0 = 0;                % Initial condition of w
% 
% % IC(2), I(2);
% vars(4).I = 0.2;                  % External input
% vars(4).u0 = -0.1;                % Inital condition of u
% vars(4).w0 = 0.1;                % Initial condition of w


%% Cycle through examples
for exampi = 1:length(vars)
    % assign variables
    u0 = vars(exampi).u0;
    w0 = vars(exampi).w0;
    I = vars(exampi).I;
    
    [u{exampi}, w{exampi}, du{exampi}, dw{exampi},...
        dU{exampi}, dW{exampi}, uNull{exampi}, wNull{exampi}] = ...
        FHNkernel(...
        epsilon,b,I,u0,w0,uvals,wvals,t);
end


%% Grouped plot
[X1, X2] = meshgrid(uvals,wvals);
U = X1(:);
W = X2(:);

us = unique(uvals);
us = linspace(min(us),max(us),100);

figure('Name','Dynamic trajectory','Position',[304 112 879 678])
subplot(2,2,1)
quiver(U,W,dU{1},dW{1})
hold on
plot(us,uNull{1},'k')
plot(us,wNull{1},'k')
plot(u{1},w{1},'.','LineWidth',2,'Color',[198 79 157]/255)
axis([min(uvals) max(uvals) min(wvals) max(wvals)])
axis square
xlabel('u')
ylabel('w')
plot(u{1}(end),w{1}(end),'o',...
    'Color','k','MarkerFaceColor',[198 79 157]/255,'MarkerSize',10)
plot(u{1}(1),w{1}(1),'o',...
    'Color',[198 79 157]/255,'MarkerFaceColor',[1 1 1],'MarkerSize',10)
mymakeaxis(gca)

subplot(2,2,2)
quiver(U,W,dU{2},dW{2})
hold on
plot(us,uNull{2},'k')
plot(us,wNull{2},'k')
plot(u{2},w{2},'LineWidth',2,'Color',[198 79 157]/255)
axis([min(uvals) max(uvals) min(wvals) max(wvals)])
axis square
xlabel('u')
ylabel('w')
plot(u{2}(end),w{2}(end),'o',...
    'Color','k','MarkerFaceColor',[198 79 157]/255,'MarkerSize',10)
plot(u{2}(1),w{2}(1),'o',...
    'Color',[198 79 157]/255,'MarkerFaceColor',[1 1 1],'MarkerSize',10)

plot(u{1}(1),w{1}(1),'o',...
    'Color',[198 79 157]/255,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',10)
mymakeaxis(gca)

subplot(2,2,3)
quiver(U,W,dU{3},dW{3})
hold on
plot(us,uNull{3},'k')
plot(us,wNull{3},'k')
plot(u{3},w{3},'.','LineWidth',2,'Color',[198 79 157]/255)
axis([min(uvals) max(uvals) min(wvals) max(wvals)])
axis square
xlabel('u')
ylabel('w')
plot(u{3}(end),w{3}(end),'o',...
    'Color','k','MarkerFaceColor',[198 79 157]/255,'MarkerSize',10)
plot(u{3}(1),w{3}(1),'o',...
    'Color',[198 79 157]/255,'MarkerFaceColor',[1 1 1],'MarkerSize',10)
mymakeaxis(gca)

subplot(2,2,4)
quiver(U,W,dU{4},dW{4})
hold on
plot(us,uNull{3},'Color',[0.7 0.7 0.7])
plot(us,uNull{4},'k')
plot(us,wNull{4},'k')
plot(u{4},w{4},'LineWidth',2,'Color',[198 79 157]/255)
axis([min(uvals) max(uvals) min(wvals) max(wvals)])
axis square
xlabel('u')
ylabel('w')
plot(u{4}(end),w{4}(end),'o',...
    'Color','k','MarkerFaceColor',[198 79 157]/255,'MarkerSize',10)
plot(u{4}(1),w{4}(1),'o',...
    'Color',[198 79 157]/255,'MarkerFaceColor',[1 1 1],'MarkerSize',10)
mymakeaxis(gca)


%% Functions
function [u, w, du, dw, dU, dW, uNull, wNull] = FHNkernel(...
    epsilon,b,I,u0,w0,uvals,wvals,t)

%% Determine vector field
[X1, X2] = meshgrid(uvals,wvals);
U = X1(:);
W = X2(:);
[dU, dW, uNull, wNull, us] = FitzHughNagumoModelVectorField2(U,W,epsilon,b,I);

%% Run simulation
[u, w, du, dw] = FitzHughNagumoModel2(t,u0,w0,epsilon,b,I,0);
    
%% Plot the trajectories, vector field and null-clines
figure('Name','Dynamic trajectory','Position',[304 112 879 394])
subplot(1,2,1)
quiver(U,W,dU,dW)
hold on
plot(us,uNull,'k')
plot(us,wNull,'k')
plot(u,w,'LineWidth',2)
axis([min(uvals) max(uvals) min(wvals) max(wvals)])
axis square
xlabel('u')
ylabel('w')
plot(u(end),w(end),'rs','MarkerFaceColor','r','MarkerSize',10)
plot(u(1),w(1),'ro','MarkerFaceColor','r','MarkerSize',10)
mymakeaxis(gca)

subplot(1,2,2)
plot(t,u,'LineWidth',2);
hold on
plot(t,w,'LineWidth',2);
axis([min(t) max(t) min(uvals) max(uvals)])
xlabel('Time')
ylabel('u/w')
mymakeaxis(gca)