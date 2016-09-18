function [mu, k, h, llike, exitflg, output] = LNP_fitter(x,y,mu_ini,k_ini,h_ini,varargin)
%% LNP_fitter
%
%   out = LNP_fitter(x,y)
%       Fits a linear-nonlinear-poisson (LNP) model to known inputs, x, and
%       recorded outputs, y, of neural data.
%
%%

%% Defaults


%% Parse inputs
Parser = inputParser;
addRequired(Parser,'x');
addRequired(Parser,'y');
addRequired(Parser,'mu_ini');
addRequired(Parser,'k_ini');
addRequired(Parser,'h_ini');
addParameter(Parser,'BasisSet',[]);
addParameter(Parser,'optMethod','ML');
addParameter(Parser,'sigma',Inf);
addParameter(Parser,'OPTIONS',[]);

parse(Parser,x,y,mu_ini,k_ini,h_ini,varargin{:})

x = Parser.Results.x;
y = Parser.Results.y;
mu_ini = Parser.Results.mu_ini;
k_ini = Parser.Results.k_ini;
h_ini = Parser.Results.h_ini;
BasisSet = Parser.Results.BasisSet;
optMethod = Parser.Results.optMethod;
sigma = Parser.Results.sigma;
OPTIONS = Parser.Results.OPTIONS;

ksz = length(k_ini);
hsz = length(h_ini);

if isempty(OPTIONS)
    OPTIONS = optimset('Display','iter');
end

if isempty(BasisSet)
    BasisSet.k = ones(size(k_ini));
    BasisSet.h = ones(size(h_ini));
end

%% Set up cross validation


%% Fit model to data
% Optimization problem
switch optMethod
    case 'ML'
        minimizant = @(p)logLike(p(1),p(2:2+ksz-1),p(2+ksz:2+ksz+hsz-1),x,y,BasisSet);
    case 'MAP'
        if isinf(sigma)
            error(['MAP fitter requires sigma input'])
        end
        minimizant = @(p)logposterior(p(1),p(2:2+ksz-1),p(2+ksz:2+ksz+hsz-1),x,y,BasisSet,sigma);
    case 'LSQ'
        minimizant = @(p)lsq(p(1),p(2:2+ksz-1),p(2+ksz:2+ksz+hsz-1),x,y,BasisSet);
end

[lparams, llike, exitflg, output] = fminsearch(minimizant,[mu_ini,k_ini,...
                                    h_ini],OPTIONS);
                      
% Unpack ML parameters
mu = lparams(1);
k = lparams(2:2+ksz-1);
h = lparams(2+ksz:2+ksz+hsz-1);

%% Functions

% Log likelihood of parameters, given data
function logL = logLike(mu,k,h,x,y,BasisSet)
% for i = 1:size(x,2)
%     Ktemp = fliplr([sum(repmat(k(i,:)',[1,size(BasisSet.k,2)]).*BasisSet.k,1) zeros(1,size(BasisSet.k,2)-1)]);
%     kX(:,i,:) = convn(x(:,i,:),Ktemp(:),'same');
% end
% kx = sum(kX,2);
% for i = 1:size(y,2)
%     Htemp = fliplr([sum(repmat(h(i,:)',[1,size(BasisSet.h,2)]).*BasisSet.h,1) zeros(1,size(BasisSet.h,2)-1)]);
% %     Htemp = fliplr([h(i,:) zeros(1,size(h,2)-1)]);
%     hY(:,i,:) = convn(y(:,i,:),Htemp(:),'same');
% end
% hy = sum(hY,2);
% ll = y .* (mu + kx + hy) - exp(mu + kx + hy);
% logL = -sum(ll(:));

for i = 1:size(x,2)
    kX(:,i,:) = filter(sum(repmat(k(i,:)',[1,size(BasisSet.k,2)]).*BasisSet.k,1),1,x(:,i,:));
end
kx = sum(kX,2);
for i = 1:size(y,2)
    hY(:,i,:) = filter([0 sum(repmat(h(i,:)',[1,size(BasisSet.h,2)]).*BasisSet.h,1)],1,y(:,i,:));
end
hy = sum(hY,2);
ll = y .* (mu + kx + hy) - exp(mu + kx + hy);
logL = -sum(ll(:));


% MAP parameters, given data
function logP = logposterior(mu,k,h,x,y,BasisSet,sigma)
% for i = 1:size(x,2)
%     Ktemp = fliplr([sum(repmat(k(i,:)',[1,size(BasisSet.k,2)]).*BasisSet.k,1) zeros(1,size(BasisSet.k,2)-1)]);
%     kX(:,i,:) = convn(x(:,i,:),Ktemp(:),'same');
% end
% kx = sum(kX,2);
% for i = 1:size(y,2)
%     Htemp = fliplr([sum(repmat(h(i,:)',[1,size(BasisSet.h,2)]).*BasisSet.h,1) zeros(1,size(BasisSet.h,2)-1)]);
% %     Htemp = fliplr([h(i,:) zeros(1,size(h,2)-1)]);
%     hY(:,i,:) = convn(y(:,i,:),Htemp(:),'same');
% end
% hy = sum(hY,2);
% ll = y .* (mu + kx + hy) - exp(mu + kx + hy);
% logprior = -[mu, k, h]*diag(sigma*ones(size([mu, k, h])))*[mu, k, h]';
% logP = -(logprior + sum(ll(:)));

for i = 1:size(x,2)
    kX(:,i,:) = filter(sum(repmat(k(i,:)',[1,size(BasisSet.k,2)]).*BasisSet.k,1),1,x(:,i,:));
end
kx = sum(kX,2);
for i = 1:size(y,2)
    hY(:,i,:) = filter(sum(repmat(h(i,:)',[1,size(BasisSet.h,2)]).*BasisSet.h,1),1,[zeros([1 1 size(y,3)]); y(1:end-1,i,:)]);
end
hy = sum(hY,2);
ll = y .* (mu + kx + hy) - exp(mu + kx + hy);
logprior = -[mu, k, h]*inv(diag(sigma*ones(size([mu, k, h]))))*[mu, k, h]';
logP = -(logprior + sum(ll(:)));

% Least Squares
function sse = lsq(mu,k,h,x,y,BasisSet)
% for i = 1:size(x,2)
%     Ktemp = fliplr([sum(repmat(k(i,:)',[1,size(BasisSet.k,2)]).*BasisSet.k,1) zeros(1,size(BasisSet.k,2)-1)]);
%     kX(:,i,:) = convn(x(:,i,:),Ktemp(:),'same');
% end
% kx = sum(kX,2);
% for i = 1:size(y,2)
%     Htemp = fliplr([sum(repmat(h(i,:)',[1,size(BasisSet.h,2)]).*BasisSet.h,1) zeros(1,size(BasisSet.h,2)-1)]);
% %     Htemp = fliplr([h(i,:) zeros(1,size(h,2)-1)]);
%     hY(:,i,:) = convn(y(:,i,:),Htemp(:),'same');
% end
% hy = sum(hY,2);
% 
% yest = exp(mu + kx + hy);
% sse = sum( (y(:) - yest(:)).^2 );

for i = 1:size(x,2)
    kX(:,i,:) = filter(sum(repmat(k(i,:)',[1,size(BasisSet.k,2)]).*BasisSet.k,1),1,x(:,i,:));
end
kx = sum(kX,2);
for i = 1:size(y,2)
    hY(:,i,:) = filter(sum(repmat(h(i,:)',[1,size(BasisSet.h,2)]).*BasisSet.h,1),1,[zeros([1 1 size(y,3)]); y(1:end-1,i,:)]);
end
hy = sum(hY,2);

yest = exp(mu + kx + hy);
sse = sum( (y(:) - yest(:)).^2 );

% 
% for j = 1:size(x,2)
%     for i = 1:size(x,1)
%         if i-length(k)+1 < 1
%             xtemp = x(1:i,j);
%             ytemp = y(1:i,j);
%             ll(i,j) = y(i,j) .* (mu + k(end-length(1:i)+1:end)*xtemp + h(end-length(1:i)+1:end)*ytemp) - exp( mu + k(end-length(1:i)+1:end)*xtemp + h(end-length(1:i)+1:end)*ytemp );
%         else
%             xtemp = x(i-length(k)+1:i,j);
%             ytemp = y(i-length(h)+1:i,j);
%             ll(i,j) = y(i,j) .* (mu + k*xtemp + h*ytemp) - exp( mu + k*xtemp + h*ytemp );
%         end
%     end
% end
% logL = -sum(ll(:));