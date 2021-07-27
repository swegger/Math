%% Prior
tau = 1000;
K = @(t1,t2)(exp(-(t1-t2).^2/tau));
t = 0:100;
[T1,T2] = meshgrid(t);
mu = zeros(length(t),1);
Sig = K(T1,T2);

%% Select a set of points
f = mvnrnd(mu,reshape(K(T1(:),T2(:)),size(T1)));


%% Select a subset of observed points and condition on those points
obsInds = true(size(t));
obsInds(21:80) = false;
tObs = t(obsInds);
fObs = f(obsInds);
tUn = t(~obsInds);

% Partition mean and covariance
muObs = mu(obsInds);
muUn = mu(~obsInds);

[T1,T2] = meshgrid(tUn);
SigUnUn = reshape(K(T1(:),T2(:)),size(T1))';
[T1,T2] = meshgrid(tUn,tObs);
SigUnObs = reshape(K(T1(:),T2(:)),size(T1))';
[T1,T2] = meshgrid(tObs,tUn);
SigObsUn = reshape(K(T1(:),T2(:)),size(T1))';
[T1,T2] = meshgrid(tObs);
SigObsObs = reshape(K(T1(:),T2(:)),size(T1))';

% Conditional mean and covaraiance
[mu_, Sig_] = conditionalGaussian(mu,Sig,fObs(:),'x_indices',obsInds);
% [mu_, Sig_] = conditionalGaussian(mu,Sig,fObs(:),'x_indices',obsInds,...
%     'SigUnUn',SigUnUn','SigObsObs',SigObsObs,'SigUnObs',SigUnObs,'SigObsUn',SigObsUn);
% mu_ = muUn + SigUnObs/SigObsObs*(fObs'-muObs);
% Sig_ = SigUnUn - SigUnObs/SigObsObs*SigObsUn;
% Sig_ = nearestSPD(Sig_);

%% Plotting
figure;
subplot(1,2,1)
imagesc(Sig_)

subplot(1,2,2)
plot(t,f,'k-')
hold on
plot(tObs,fObs,'ro')
plot(tUn,mu_,'bx')
% plot(tUn,mvnrnd(mu_,Sig_),'-.')
% plot(tUn,mvnrnd(mu_,Sig_),'-.')
