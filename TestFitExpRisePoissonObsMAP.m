%% TestFitExpRisePoissonObsMAP
%
%   Tests the fitter on known fake data
%
%%

%% Variables
trials = 100;
xs = 0:50:500;
theta = [10,15,100,200];
f = @(x,theta)(pieceWiseExpRise(x,theta));

%% Generate data
x = repmat(xs,[trials,1]);
lambda = f(x,theta);
y = poissrnd(lambda);

%% Fit the data
theta0 = [10*randn,20*randn,200*rand,300*rand];
while theta0(1)+theta0(2) < 0
    theta0(1:2) = [10*randn,20*randn];
end
[thetaFit, logPosterior] = FitExpRisePoissonObsMAP(x(:),y(:),theta0,'mu',NaN,'sig',NaN);
% [thetaFitG, logPosteriorG] = FitGaussPoissonObsMAP(x(:),y(:),[2,1000,200],'mu',NaN,'sig',NaN);
[thetaFitc, logPosteriorc] = FitConstantPoissonObsMAP(x(:),y(:));%,'theta_init',[100],'a',NaN,'b',NaN);

%% Plot the results
figure('Name','Fit results')
plot(xs,f(xs,theta),'k');
hold on
plot(xs,f(xs,thetaFit),'r');
plot(xs,thetaFitc*ones(size(xs)),'r--');
% for i = 1:length(xs)
%     h(4) = errorbar(xs(i),mean(y(x == xs(i))),std(y(x==xs(i)))/sqrt(sum(x==xs(i))),'ko');
% end
plot(x,y,'.','Color',[0.6 0.6 0.6]);
xlabel('Input')
ylabel('Output')
% legend(h,{'Actual tuning','Fit tuning','Fit constant','Mean +/- ste','Data'})
mymakeaxis(gca)

% figure('Name','Log likelihood of models, given data')
% bar([1 2 3],[logPosterior, logPosteriorG, logPosteriorc],'k')