%% TestFitGaussPoissonObsMAP
%
%   Tests the fitter on known fake data
%
%%

%% Variables
trials = 1000;
xs = 600:100:1000;
theta = [1,1200,600];
mu = [1;1200;200];
sig = [Inf,0,0,;0,Inf,0;0,0,Inf];
f = @(x,theta)(theta(1)*exp(-(x-theta(2)).^2./(2*theta(3)^2)));

%% Generate data
x = xs(ceil(length(xs)*rand(trials,1)));
lambda = f(x,theta);
y = poissrnd(lambda);

%% Fit the data
[thetaFit, logPosterior] = FitGaussPoissonObsMAP(x(:),y(:),[2,1000,200],'mu',NaN,'sig',NaN);
[thetaFitL, logPosteriorL] = FitLinearPoissonObsMAP(x(:),y(:),[100,10],'mu',NaN,'sig',NaN);
[thetaFitc, logPosteriorc] = FitConstantPoissonObsMAP(x(:),y(:));%,2,'mu',NaN,'sig',NaN);

%% Plot the results
figure('Name','Fit results')
h(1) = plot(600:1000,f(600:1000,theta),'k');
hold on
h(2) = plot(600:1000,f(600:1000,thetaFit),'r');
h(3) = plot(600:1000,thetaFitc*ones(size(600:1000)),'r--');
for i = 1:length(xs)
    h(4) = errorbar(xs(i),mean(y(x == xs(i))),std(y(x==xs(i)))/sqrt(sum(x==xs(i))),'ko');
end
h(5) = plot(x,y,'.','Color',[0.6 0.6 0.6]);
xlabel('Input')
ylabel('Output')
legend(h,{'Actual tuning','Fit tuning','Fit constant','Mean +/- ste','Data'})
mymakeaxis(gca)

figure('Name','Log likelihood of models, given data')
bar([1 2 3],[logPosterior, logPosteriorL, logPosteriorc],'k')