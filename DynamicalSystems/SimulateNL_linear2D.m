%% FitzHughNagumoToyExample
%
%   Toy example of the FitzHugh-Nagumo 2D dynamical model to determine if
%   observed trajectories can be explained by a simple model.
%
%%

% Variables
dt = 0.001;              % Time step of simulation
A = [-1 0; 0 -1];
gain = 0.1;
F{1} = @(x)( x(:,1) ./ sqrt(1+x(:,1).^2) );
F{2} = @(x)( gain*exp(x(:,2)) );
I = 0;                  % External input
B = [0 0; 0 0];
u0 = [1; 1];                % Inital condition of u
uvals = -3:0.1:3;
kick = [0 -2];

tkick = 4*[0.6 0.7 0.8 0.9 1];

colors = [     0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880];%...
%     0.3010    0.7450    0.9330;...
%     0.6350    0.0780    0.1840];
colors = flipud(colors);

%% Determine vector field
[X1, X2] = meshgrid(uvals);
U = [X1(:), X2(:)];
[dX, XNull] = NL_linear2DVectorField(U,A,F,I);

%% Run simulation and "kick" the partical to a new position at different
% times
for i = 1:length(tkick)
    % Run the model up until the kick
    t{i} = 0:dt:tkick(i);
    [u{i}, du{i}] = NL_linear2D(t{i},u0,A,F,I,B);
    
    % Run the model following the kick
    kicktemp = -kick;%*w{i}(end);
    [upost{i}, dupost{i}] = NL_linear2D(t{i},u{i}(end)+kicktemp(1),A,F,I,B);
end

%% Measure the "speed" of activity changes along u and w
for i = 1:length(tkick)
    q{i} = (du{i}.^2 + dw{i}.^2)/2;
    qpost{i} = (dupost{i}.^2 + dwpost{i}.^2)/2;
end

%% Plot the trajectories, vector field and null-clines
figure('Name','Dynamic trajectory','Position',[250 247 560 420])
quiver(U,W,dU,dW)
hold on
plot(uvals,uNull,'k')
plot(uvals,wNull,'k')
for i = 1:length(tkick)
    plot(u{end-i+1},w{end-i+1},'LineWidth',2,'Color',colors(i,:))
    plot(upost{end-i+1},wpost{end-i+1},'LineWidth',2,'Color',colors(i,:))
end
axis([min(uvals) max(uvals) min(wvals) max(wvals)])
xlabel('u')
ylabel('w')

figure('Name','Dynamics of u and w, pre-kick','Position',[812 246 560 420])
for i = 1:length(tkick);
    subplot(2,1,1)
    plot(t{length(tkick)-i+1},u{length(tkick)-i+1},'LineWidth',2,'Color',colors(i,:))
    hold on
    ylabel('u')
    subplot(2,1,2)
    plot(t{length(tkick)-i+1},w{length(tkick)-i+1},'LineWidth',2,'Color',colors(i,:))
    hold on
    ylabel('w')
end
xlabel('time')

figure('Name','Dynamics of u and w, post-kick','Position',[1372 246 560 420])
for i = 1:length(tkick);
    subplot(2,1,1)
    plot(t{length(tkick)-i+1},upost{length(tkick)-i+1},'LineWidth',2,'Color',colors(i,:))
    hold on
    ylabel('u')
    subplot(2,1,2)
    plot(t{length(tkick)-i+1},wpost{length(tkick)-i+1},'LineWidth',2,'Color',colors(i,:))
    hold on
    ylabel('w')
end
xlabel('time')

figure('Name','Speed profile')
for i = 1:length(tkick)
    subplot(2,1,1)
    plot(t{length(tkick)-i+1},q{length(tkick)-i+1},'LineWidth',2,'Color',colors(i,:));
    hold on
    subplot(2,1,2)
    plot(t{length(tkick)-i+1},qpost{length(tkick)-i+1},'LineWidth',2,'Color',colors(i,:));
    hold on
end
subplot(2,1,1)
ylabel('q(t), pre-kick')
subplot(2,1,2)
xlabel('time')
ylabel('q(t), post-kick')
