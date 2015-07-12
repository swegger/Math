%% TestFitMixtrureModelEM
%
%   Tests the fitter
%
%
%%

% Parameters
N = 30000;
D = 2;
tolerance = 0.01;
K = 3;
w = [0.1 0.4 0.5];
support = [-20 20];
noiseModel.type = 'uniform';
nBackgroundPoints = 10;

switch D
    case 2
        % 3D mixture model
        mu = [-5 -5; 5 5; -5 5];
        sig(1,1,1) = 0.5;
        sig(1,2,1) = 0.1;
        sig(2,1,1) = 0.1;
        sig(2,2,1) = 0.5;
        sig(1,1,2) = 0.5;
        sig(1,2,2) = 0.1;
        sig(2,1,2) = 0.1;
        sig(2,2,2) = 0.5;
        sig(1,1,3) = 0.5;
        sig(1,2,3) = 0.1;
        sig(2,1,3) = 0.1;
        sig(2,2,3) = 0.5;
        x = nan(N,D);
        dx = 0.1;
        xvec = min(mu)-5*abs(max(reshape(sig,numel(sig),1))*min(mu)):dx:max(mu)+5*abs(max(reshape(sig,numel(sig),1))*max(mu));
        
        w_init = [0.5 0.5 0.5];
        mu_init = [-1 -1; 0 0; -1 1];
        sig_init(1,1,1) = 1;
        sig_init(1,2,1) = 0;
        sig_init(2,1,1) = 0;
        sig_init(2,2,1) = 1;
        sig_init(1,1,2) = 1;
        sig_init(1,2,2) = 0;
        sig_init(2,1,2) = 0;
        sig_init(2,2,2) = 1;
        sig_init(1,1,3) = 1;
        sig_init(1,2,3) = 0;
        sig_init(2,1,3) = 0;
        sig_init(2,2,3) = 1;
        
    case 3
        % 3D mixture model
        mu = [-5 -5 -5; 5 5 5; -5 5 5];
        sig(1,1,1) = 0.5;
        sig(1,2,1) = 0.1;
        sig(1,3,1) = 0.1;
        sig(2,1,1) = 0.1;
        sig(2,2,1) = 0.5;
        sig(2,3,1) = 0.1;
        sig(3,1,1) = 0.1;
        sig(3,2,1) = 0.1;
        sig(3,3,1) = 0.5;
        sig(1,1,2) = 0.5;
        sig(1,2,2) = 0.1;
        sig(1,3,2) = -0.1;
        sig(2,1,2) = 0.1;
        sig(2,2,2) = 0.5;
        sig(2,3,2) = 0.1;
        sig(3,1,2) = -0.1;
        sig(3,2,2) = 0.1;
        sig(3,3,2) = 0.5;
        sig(1,1,3) = 0.5;
        sig(1,2,3) = 0.1;
        sig(1,3,3) = 0.1;
        sig(2,1,3) = 0.1;
        sig(2,2,3) = 0.5;
        sig(2,3,3) = -0.1;
        sig(3,1,3) = 0.1;
        sig(3,2,3) = -0.1;
        sig(3,3,3) = 0.5;
        x = nan(N,D);
        dx = 0.1;
        xvec = min(mu)-5*abs(max(reshape(sig,numel(sig),1))*min(mu)):dx:max(mu)+5*abs(max(reshape(sig,numel(sig),1))*max(mu));
        
        w_init = [0.5 0.5 0.5];
        mu_init = [-1 0 0; 1 0 0; 0 0 0];
        sig_init(1,1,1) = 1;
        sig_init(1,2,1) = 0;
        sig_init(1,3,1) = 0;
        sig_init(2,1,1) = 0;
        sig_init(2,2,1) = 1;
        sig_init(2,3,1) = 0;
        sig_init(3,1,1) = 0;
        sig_init(3,2,1) = 0;
        sig_init(3,3,1) = 1;
        sig_init(1,1,2) = 1;
        sig_init(1,2,2) = 0;
        sig_init(1,3,2) = 0;
        sig_init(2,1,2) = 0;
        sig_init(2,2,2) = 1;
        sig_init(2,3,2) = 0;
        sig_init(3,1,2) = 0;
        sig_init(3,2,2) = 0;
        sig_init(3,3,2) = 1;
        sig_init(1,1,3) = 1;
        sig_init(1,2,3) = 0;
        sig_init(1,3,3) = 0;
        sig_init(2,1,3) = 0;
        sig_init(2,2,3) = 1;
        sig_init(2,3,3) = 0;
        sig_init(3,1,3) = 0;
        sig_init(3,2,3) = 0;
        sig_init(3,3,3) = 1;

end

% Throw coin to see which Gaussian to use
for i = 1:K
    rn = rand(N,1);
    if i == 1
        z(:,i) = rn <= w(1);
    elseif i == K
        z(:,i) = rn > sum(w(1:end-1));
    else
        z(:,i) = rn > sum(w(1:i-1)) & rn <= sum(w(1:i));
    end
end

% Generate points
for i = 1:K
    x(z(:,i),:) = mvnrnd(mu(i,:),sig(:,:,i),sum(z(:,i)));   %(sig(:,:,i)*randn(D,sum(z(:,i))) + repmat(mu(i,:)',1,sum(z(:,i))))';
end
x = x(~isnan(x(:,1)),:);
x = [x; (support(2)-support(1))*rand(nBackgroundPoints,D) + support(1)];
%x(isnan(x(:,1)),:) = (support(2)-support(1))*rand(sum(isnan(x(:,1))),D) + support(1);        % Replace NaNs with a random point sampled from a uniform distribution with volume (support(2)-support(1))^D

% Fit the model to the data
[w_ml, mu_ml, sig_ml, logl, gamma] = FitMixtureModelEM(x,K,w_init,mu_init,sig_init,tolerance,'Batch',1000,'percentFit',50,'noiseModel',noiseModel);

% Plot the data and the fit
switch D
    case 1
        figure
        n = hist(x,xvec);
        pdata = n/sum(n)/dx;
        pfit = GaussMix(w_ml,mu_ml,sig_ml,xvec(:));
        
        figure
        bar(xvec,pdata,'k')
        hold on
        plot(xvec,pfit,'r','LineWidth',2)
        
    case 2
        figure
        scatter(x(:,1),x(:,2),10,gamma);
        hold on
        for i = 1:K
            plot(mu(i,1),mu(i,2),'k*')
            [v, d] = eig(sig(:,:,i));
            [ex, ey] = generalEllipse(3*sqrt(d(end,end)),3*sqrt(d(end-1,end-1)),'Center',mu(i,:),'theta',0:pi/180:2*pi,'phi',atan(v(2,2)/v(1,2)));
            plot(ex,ey,'k','LineWidth',2)
            plot(mu_ml(i,1),mu_ml(i,2),'m*')
            [v, d] = eig(sig_ml(:,:,i));
            [ex, ey] = generalEllipse(3*sqrt(d(end,end)),3*sqrt(d(end-1,end-1)),'Center',mu_ml(i,:),'theta',0:pi/180:2*pi,'phi',atan(v(2,2)/v(1,2)));
            plot(ex,ey,'m','LineWidth',2)
        end
        legend('data','True means','True covariance','Fit means','Fit covariance')
end