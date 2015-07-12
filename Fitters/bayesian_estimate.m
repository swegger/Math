function estimate = bayesian_estimate( ...
    domain, ...
    measurements, ...
    weberMeasure, ...
    support)
%% Using Integral functions
%Some of the inputs here aren't used, but can be.
wm = weberMeasure;

alpha = 1./integral(@(s)posteriorNum(s,measurements,wm), support(1), support(2), ...
    'ArrayValued',true);
estimate = integral(@(s)alpha.*s.*posteriorNum(s,measurements,wm), support(1), support(2), ...
    'ArrayValued',true);

% posteriorBeta = arrayfun(@(s) alpha.*posteriorNum(s,measurements,wm),domain, ...
%     'UniformOutput',false);



end

function numValue = posteriorNum(s,m,wm)
%% Using multiprod(), assuming a single value of s passed in at a time.
%This is the simplest version of the calculation. Does not allow the
%measurements to covary.
sigmaDet = ((s*wm).^2)^size(m,2);
numValue = 1./sqrt(2*pi*sigmaDet)*exp(-.5*sum((s-m).^2/(s*wm)^2,2));
end

function numValue = posteriorNum2(s,m,wm)
%% Using multiprod(), assuming a single value of s passed in at a time.
%Generalized multivariate Gaussian, can handle covarying measurements.
sigmaMat = diag(repmat((wm*s).^2,size(m,2),1));

exponExpres = multiprod(s-m,inv(sigmaMat));
exponExpres = multiprod(exponExpres,s-m,2,2);
numValue = 1./sqrt(2*pi*det(sigmaMat))*exp(-0.5*exponExpres);
end

function numValue = posteriorNum3(s,m,wm)
%% Using mmx(), assuming a single value of s passed in at a time.
%Generalized multivariate Gaussian, can handle covarying measurements.
m = reshape(m,[size(m,2) 1 size(m,1)]);
sigmaMat = diag((wm*s).^2*ones(size(m,1),1));
invSigmaMat = inv(sigmaMat);
detSigmaMat = det(sigmaMat);
s = repmat(s,size(m,1),1,size(m,3));

exponExpres = mmx('mult',permute(s-m,[2 1 3 4]),invSigmaMat);
exponExpres = mmx('mult',exponExpres,s-m);
exResult = exp(-0.5*exponExpres);

numValue = 1./sqrt(2*pi*detSigmaMat).*exResult;
numValue = squeeze(reshape(numValue,[size(numValue,3) size(numValue,4) 1 1]));
end

function numValue = posteriorNum4(s,m,wm)
%% This function can generate output over arrays of s and m independently.
%Generalized multivariate Gaussian, can handle covarying measurements.
%Unfortunately, the integrate() function can't make use of it. Seems to be
%a bit slower without this function being utilized. Could make my own
%integration function...
m = reshape(m,[size(m,2) 1 size(m,1) 1]);
m = repmat(m,1,1,1,numel(s));

sigmaMat = repmat(diag(ones(size(m,1),1)),[1 1 1 numel(s)]);
sigmaMat = sigmaMat.*repmat(reshape((wm*s).^2,[1 1 1 numel(s)]),size(m,2),size(m,2));

invSigmaMat = zeros(size(sigmaMat));
for i = 1:size(sigmaMat,4)
    invSigmaMat(:,:,1,i) = inv(sigmaMat(:,:,1,i));
end

detSigmaMat = zeros(1,1,size(m,3),size(m,4));
for i = 1:numel(s)
    detSigmaMat(1,1,:,i) = det(sigmaMat(:,:,1,i));
end

s = reshape(s,[1 1 1 numel(s)]);
s = repmat(s,size(m,1),1,size(m,3),1);

exponExpres = mmx('mult',permute(s-m,[2 1 3 4]),invSigmaMat);
exponExpres = mmx('mult',exponExpres,s-m);
exResult = exp(-0.5*exponExpres);

numValue = 1./sqrt(2*pi*detSigmaMat).*exResult;
numValue = squeeze(reshape(numValue,[size(numValue,3) size(numValue,4) 1 1]));
end

%% Doing integration "by hand"
% %% Calculate likelihood functions for each sample
% posteriorStart = find(domain >= support(1),1,'first');
% posteriorEnd = find(domain <= support(2),1,'last');
% numelSupport = posteriorEnd - posteriorStart + 1;
% effDomain = domain(posteriorStart:posteriorEnd);
% effSigmaM = weberMeasure*effDomain;
% likelihood = zeros(numel(domain),n_simulations);
% likelihood(posteriorStart:posteriorEnd,:) = ones(numelSupport,n_simulations);
% %%
% for i = 1:size(measurements,2)
% likelihood(posteriorStart:posteriorEnd,:) = ...
%     likelihood(posteriorStart:posteriorEnd,:).*...
%     exp(-(repmat(effDomain,1,n_simulations)- ...
%     repmat(measurements(:,i)',numelSupport,1)).^2./ ...
%     (2*repmat(effSigmaM.^2,1,n_simulations)))./ ...
%     sqrt(2*pi*repmat(effSigmaM.^2,1,n_simulations));
% end
% 
% %WARNING: if you give this function a measurement_sigma that is much
% %smaller than the variance of the distribution used to generate the
% %measurements, the the likelihood may be equal to zero everywhere where the
% %prior is nonzero. Then the likelihood*prior expression can return a vector
% %of zeros, and the estimate will be a NaN.
% 
% %Generate posterior distribution for each likelihood function
% posterior = likelihood;
% posterior = posterior./repmat(sum(posterior),length(domain),1);
% 
% %make mmse estimate = E[X|Y = y];
% estimate = expectation(posterior,domain);