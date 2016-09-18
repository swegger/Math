%% TestFitExponential
%
%
%%

% Parameters
offset = 0;
tau = 0.01;
x = 0:1000;
y = offset + 0.1*exp(tau*x);

%% Fit the data
[theta, feval] = FitExponential(x,y,[0.0001 0.1]);

%% Plot the results
figure
plot(x,y,'k')
hold on
plot(x,theta(2) + exp(theta(1)*x),'r');