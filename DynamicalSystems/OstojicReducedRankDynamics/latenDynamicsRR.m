function [dkappas, kappas, sigmaTildes, deltas, vs] = latenDynamicsRR(tau,us,overlaps,kappas0,sigmas,varargin)
%% latentDynamicsRR
%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'tau')
addRequired(Parser,'us')
addRequired(Parser,'overlaps')
addRequired(Parser,'kappas0')
addRequired(Parser,'sigmas')
addParameter(Parser,'eps',1e-10)
addParameter(Parser,'nsteps',1000)

parse(Parser,tau,us,overlaps,kappas0,sigmas,varargin{:})

tau = Parser.Results.tau;
us = Parser.Results.us;
overlaps = Parser.Results.overlaps;
kappas0 = Parser.Results.kappas0;
sigmas = Parser.Results.sigmas;
eps = Parser.Results.eps;
nsteps = Parser.Results.nsteps;

%% Initialize
N = length(kappas0);
M = size(us,1);
T = size(us,2);

Orec = overlaps(:,1:N);
Oinput = overlaps(:,N+1:end);
if size(Oinput,2) ~= M
    error('overlaps should be N by N+M')
end

dvs = nan(M,T);
vs = nan(M,T);

deltas = nan(1,T);
sigmaTildes = nan(N,N+M);

dkappas = nan(N,T);
kappas = nan(N,T);

%% Run simulation

vs(:,1) = us(:,1);
kappas(:,1) = kappas0;
for ti = 2:T
    dvs(:,ti) = -vs(:,ti-1) + us(:,ti);
    vs(:,ti) = vs(:,ti-1) + dvs(:,ti)/tau;
    
    deltas(ti) = computeDelta(kappas(:,ti-1),us(:,ti),sigmas);
    sigmaTildes(:,:,ti) = overlaps .* deltaMean(deltas(ti),eps,nsteps);
    
    dkappas(:,ti) = -kappas(:,ti-1) + sigmaTildes(:,1:N,ti)*kappas(:,ti-1) + ...
        sigmaTildes(:,N+1:end,ti)*vs(:,ti);
    kappas(:,ti) = kappas(:,ti-1) + dkappas(:,ti)/tau;
end


%% Functions

%% computeDelta
function delta = computeDelta(kappas,us,sigmas)
    %%
    a = sum(kappas.^2 .* sigmas(1:size(kappas,1)).^2) + ...
        sum(us.^2 .* sigmas(size(kappas,1)+1:size(kappas,1)+size(us,1)).^2);
    delta = sqrt(a);
    
%% deltaMean
function mu = deltaMean(deltas,eps,nsteps)
    %%
    if deltas>0
        zs = linspace(atanh(-1+eps)/deltas,atanh(1-eps)/deltas,nsteps);
        dz = zs(2)-zs(1);
        mu = (1/sqrt(2*pi)) * trapz( dz * exp(-zs.^2/2) .* ...
            (1 - tanh(deltas*zs).^2) );
    else
        mu = 0;
    end    
    if any(mu(:) > 1)
        error('Check here')
    end