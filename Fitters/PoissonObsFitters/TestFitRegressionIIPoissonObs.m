%% TestFitLinearPoissonObsMAP
%
%   Tests the fitter on known fake data
%
%%

%% Variables
trials = 100;
theta = [0.5, 0.5, 0];

%% Generate data
x = poissrnd(theta(1),[trials,1]);
lambda = x*theta(2) + theta(3);
lambda(lambda < 0) = 0;
y = poissrnd(lambda);

%% Fit the data
[thetaFit, logl] = FitRegressionIIPoissonObs(x(:),y(:),[2 1 3]);

%% Plot the results
figure('Name','Fit results')
h(1) = plot(x,y,'k.');
hold on
xlabel('Input')
ylabel('Output')
mymakeaxis(gca)