function SimShuntingInhibition(t,N,varargin)
%%
%
%
%
%%

%% Defaults
coherence_default = randsample([20 60 100]/100,100,true);

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'t')
addRequired(Parser,'N')
addParameter(Parser,'tau',100)
addParameter(Parser,'coherence',coherence_default)
addParameter(Parser,'A',NaN)
addParameter(Parser,'B',NaN)
addParameter(Parser,'u',NaN)
addParameter(Parser,'Sigma',NaN)
addParameter(Parser,'initSigma',NaN)
addParameter(Parser,'baseline',log(10))

parse(Parser,t,N,varargin{:})

t = Parser.Results.t;
N = Parser.Results.N;
tau = Parser.Results.tau;
coherence = Parser.Results.coherence;
A = Parser.Results.A;
B = Parser.Results.B;
u = Parser.Results.u;
Sigma = Parser.Results.Sigma;
initSigma = Parser.Results.initSigma;
baseline = Parser.Results.baseline;

if any(isnan(A(:)))
    A = ones(N,1);
end

if any(isnan(B(:)))
    B = randn(N,1);
end

if any(isnan(u(:)))
    uCoh = setUdefault(t,coherence,N);
    u = zeros(N,length(t),length(coherence));
    u(:,t >= 0 & t <= 600,:) = 1;
end

if any(isnan(Sigma))
    Sigma = diag(ones(N,1))/10;
end

if any(isnan(initSigma))
    initSigma = diag(ones(N,1))/10;
end

%% Run simulation
v = nan([N,length(t),length(coherence)]);
dv = nan([N,length(t),length(coherence)]);

v(:,1,:) = initSigma*randn(N,length(coherence));
dv(:,1,:) = ( -v(:,1,:) + repmat(B,[1,1,size(u,3)]).*u(:,1,:) + ...
        permute(Sigma*randn(N,length(coherence)),[1,3,2]) )/tau;

for ti = 2:length(t)
    dv(:,ti,:) = ( -repmat(A,[1,1,size(u,3)]).*uCoh(:,ti,:).*v(:,ti-1,:) + repmat(B,[1,1,size(u,3)]).*u(:,ti,:) + ...
        permute(Sigma*randn(N,length(coherence)),[1,3,2]) )/tau;
    
    v(:,ti,:) = v(:,ti-1,:) + dv(:,ti,:);
end
% counts = double(exp(v) > 0.9);
counts = poissrnd(exp(v + baseline));
% counts(counts > 1) = 1;


%% Find mean and covariance of counts for each neuron
cohs = unique(coherence);
for ni = 1:N
    for ci = 1:length(cohs)
        m(:,ni,ci) = mean(counts(ni,:,coherence == cohs(ci)),3);
        va(:,ni,ci) = var(counts(ni,:,coherence == cohs(ci)),[],3);
        C(:,:,ni,ci) = cov(permute(counts(ni,:,coherence == cohs(ci)),[3,2,1]));
    end
end


%% Fit to data


%% Plotting
figure('Name','Lambda over time')

%% Functions

%% setUdefault
function u = setUdefault(t,coherence,N)
    %%
    u = ones(N,length(t),length(coherence));
    utemp = repmat(permute(coherence,[1,3,2]),[N,sum(t>=0 & t<=600),1]);
    u(:,t >= 0 & t <= 600,:) = 1 - 0.7*utemp;
    
%%