 function [w, mu, sig, logl, gamma] = FitMixtureModelEM(x,k,w_init,mu_init,sig_init,tolerance,varargin)
%% FitMixtureModelEM
%
%   [w, mu, sig, logl, gamma] = FitMixtureModelEM(x,k,mu_init,sig_init,tolerance)
%
%   Fits a mixture of k Gaussians to the data set x using the
%   expectation-maximumization (EM) algorithm and starting from the
%   parameters w_init, mu_init, and sig_int for the weights, means and
%   covariance, respectively. Iteration ends once the change in likelihoood
%   is less than tolerance. Returns, the maximum likelihood parameters w,
%   mu and sigma, their log likelihood given the data (logl), the posterior
%   probability that each point belongs to each Gaussian given
%   the data, p(z=k|x).
%
%   [w, mu, sig, logl, gamma] = FitMixtureModelEM(...,'Batch',batchsz)
%
%   Performs the calculation of likelihoods/probabilities of the Gaussian
%   mixture in batches of batchsz.
%
%   [w, mu, sig, logl, gamma] = FitMixtureModelEM(...,'percentFit',percent)
%
%   Fits the model to percent% of the data, randomly chosen without
%   replacement.
%
%
%%

% Defaults
noiseModel_default.type = 'none';

% Parse inputs
Parser = inputParser;
addRequired(Parser,'x');
addRequired(Parser,'k');
addRequired(Parser,'w_init');
addRequired(Parser,'mu_init');
addRequired(Parser,'sig_init');
addRequired(Parser,'tolerance');
addParameter(Parser,'Batch',0);
addParameter(Parser,'percentFit',100);
addParameter(Parser,'verbose','logLikelihood');
addParameter(Parser,'noiseModel',noiseModel_default);

parse(Parser,x,k,w_init,mu_init,sig_init,tolerance,varargin{:})

x = Parser.Results.x;
k = Parser.Results.k;
w_init = Parser.Results.w_init;
mu_init = Parser.Results.mu_init;
sig_init = Parser.Results.sig_init;
tolerance = Parser.Results.tolerance;
Batch = Parser.Results.Batch;
percentFit = Parser.Results.percentFit;
verbose = Parser.Results.verbose;
noiseModel = Parser.Results.noiseModel;

% Initialize the weights, means and covariance matrices of the fit
w = w_init;
mu = mu_init;
sig = sig_init;
logl = 0;
delta = 10000000;

% Extract noise model
switch noiseModel.type
    case 'none'
        v = 0;
    
    case 'uniform'
        v = 1./prod(max(x,[],1)-min(x,[],1));       % Uniform background noise over a volume covering the range of x values
        
    case 'custom'
        v = noiseModel.model(x,noiseModel.params);      % Custom noise model specified the the user-defined function noiseModel.model and parameters noiseModel.params
        
end

if percentFit/100 < 1
    xtemp = datasample(x,ceil(percentFit/100*size(x,1)),1,'Replace',false);
else
    xtemp = x;
end

N = size(xtemp,1);
D = size(xtemp,2);

% Iterate until change in log likelihood is below tolerance
while delta >= tolerance    
    
    % Set old parameter values
    loglold = logl;
    wold = w;
    muold = mu;
    sigold = sig;
    
    % Expectation step
    gamma = responsibilities(xtemp,wold,muold,sigold,Batch,v);
%     gamma(end,~any(gamma)) = realmin('double');
%     gamma = gamma./repmat(sum(gamma,2),1,size(gamma,2));
    
    % Maximization step
    Nk = nansum(gamma,1);
    sigtemp = zeros(size(sig));
    for i = 1:k
        if sum(gamma(:,i)) == 0
            mu(i,:) = muold(i,:);
            sig(:,:,i) = sigold(:,:,i);
        else
            mu(i,:) = 1/Nk(i)*sum(repmat(gamma(:,i),[1 D]).*xtemp,1);
            for j = 1:N
                sigtemp(:,:,i) = sigtemp(:,:,i) + gamma(j,i)*(xtemp(j,:)-mu(i,:))'*(xtemp(j,:)-mu(i,:));
            end
            sig(:,:,i) = 1/Nk(i)*sigtemp(:,:,i);
            [~, p] = chol(sig(:,:,i));
            if p ~= 0
                disp('Warning: EM sigma is not SPD; finding nearest')
                sig(:,:,i) = nearestSPD(sig(:,:,i));
            end
            
            zeroind = find(diag(sig(:,:,i)) == 0);
            sig(zeroind,zeroind,i) = realmin('double');
            
            w(i) = Nk(i)/N;
        end
    end
    
    % Compute log likelihood
    Like = GaussMix(w,mu,sig,xtemp,'Batch',Batch);
    Like(Like == 0) = realmin('double');
    logl = sum(log(Like));
    switch verbose
        case 'logLikelihood'
            disp(logl)
    end
    
    % Change in loglikelihood
    delta = abs(logl - loglold);
    
end

% Calculate the gamma values for each all data points using ML parameters
gamma = responsibilities(x,w,mu,sig,Batch,v);


%% Functions
% Responsibilites
function gamma = responsibilities(x,w,mu,sig,Batch,v)

mix = GaussMix(w,mu,sig,x,'Batch',Batch) + v;
mix(mix == 0) = realmin('double');
for i = 1:length(w)
    gamma(:,i) = w(i)*mvnpdf(x,mu(i,:),sig(:,:,i))./mix;
end