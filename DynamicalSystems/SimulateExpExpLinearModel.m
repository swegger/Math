%% SimulateExpExpLinearModel
%
%   Toy example of the FitzHugh-Nagumo 2D dynamical model to determine if
%   observed trajectories can be explained by a simple model.
%
%%

% Variables
dt = 0.001;              % Time step of simulation
Ainit = [-1 0; 0 -6];         %
B = [0 0; 0 0];
uvals = -3:0.1:3;
wvals = -3:0.1:3;

ts = [0.6 0.7 0.8 0.9 1];

%load /Volumes/jazlab-5/Projects/RSNGcc_Trigger/B/FitResults/B_FitAllOutput_EffectorType0_20151003.mat wm wp
load ~/Projects/RSNGcc_Trigger/B/FitResults/B_FitAllOutput_EffectorType0_20151003.mat wm wp
TAexpectation.On = true;
TAexpectation.method = 'numerical';
TAexpectation.trialtypes = [2];
TAexpectation.ts_vec = ts*1000;
TAexpectation.simtrials = 10000;
TAexpectation.Support = [600 1000];
TAexpectation.wm = wm(2);
TAexpectation.wp = wp(2);
method_opts.dx = 0.01;
[ta, ta_std] = ta_expectation3(TAexpectation.ts_vec,wm(2),1,10,'Type','BLS','method_options',method_opts.dx,'method',TAexpectation.method,'trials',TAexpectation.simtrials,'wp',0,'Support',TAexpectation.Support);

for i = 1:length(ts)
    A{i} = Ainit*max(ta)/ta(i);
%     A{i}(2,2) = Ainit(2,2)*(max(ta)/ta(i));
%     A{i}(1,1) = Ainit(1,1)*(max(ta)/ta(i));
end

u0 = 0.5 + log(ts/max(ts));                % Inital condition of u
w0 = 1.5 + log(ts/max(ts));                % Initial condition of w

inputmu = 0.200;
inputsig = 0.07;
inputamp = 10;

colors = [     0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880];%...
%     0.3010    0.7450    0.9330;...
%     0.6350    0.0780    0.1840];
colors = flipud(colors);


%% Create temporal input pattern
ttemp = 0:dt:max(ts);
input = inputamp * exp( -(ttemp - inputmu).^2 ./ (2*inputsig^2) );
I = repmat(input,2,1);
I(2,:) = 2*I(2,:);

%% Determine vector field
[X1, X2] = meshgrid(uvals,wvals);
U = X1(:);
W = X2(:);
[dU, dW, uNull, wNull] = ExpExpLinearModelVectorField(U,W,A{end},I);

%% Run simulation and "kick" the partical to a new position at different times
for i = 1:length(ts)
    % Run the model up until the kick
    t{i} = 0:dt:ts(i);
    Itemp = I*sqrt(max(ta)/ta(i));
    [u{i}, w{i}, du{i}, dw{i}] = ExpExpLinearModel(t{i},u0(i),w0(i),A{i},B,Itemp);
    
end

%% Measure the "speed" of activity changes along u and w
for i = 1:length(ts)
    q{i} = (du{i}.^2 + dw{i}.^2)/2;
end

%% Plot the trajectories, vector field and null-clines
figure('Name','Input')
for i = 1:length(ts)
    plot(t{i},I(1,1:length(t{i})),'Color',colors(i,:),'LineWidth',2)
    hold on
end
xlabel('Time (s)')
ylabel('Input');
mymakeaxis(gca)

figure('Name','Dynamic trajectory','Position',[250 247 560 420])
quiver(U,W,dU,dW)
hold on
plot(uvals,uNull,'k')
plot(uvals,wNull,'k')
for i = 1:length(ts)
    plot(u{end-i+1},w{end-i+1},'LineWidth',2,'Color',colors(i,:))
end
axis square
axis([min(uvals) max(uvals) min(wvals) max(wvals)])
xlabel('u')
ylabel('w')
mymakeaxis(gca);

figure('Name','Dynamics of u and w','Position',[812 246 560 420])
for i = 1:length(ts);
    subplot(2,1,1)
    plot(t{length(ts)-i+1},u{length(ts)-i+1},'LineWidth',2,'Color',colors(i,:))
    hold on
    ylabel('u')
    subplot(2,1,2)
    plot(t{length(ts)-i+1},w{length(ts)-i+1},'LineWidth',2,'Color',colors(i,:))
    hold on
    ylabel('w')
end
xlabel('time')

figure('Name','Speed profile')
for i = 1:length(ts)
    plot(t{length(ts)-i+1},q{length(ts)-i+1},'LineWidth',2,'Color',colors(i,:));
    hold on
end
ylabel('q(t)')
xlabel('time')
mymakeaxis(gca)