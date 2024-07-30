function [mu_, Sig_, muUn, SigUnObs, SigObsObs, SigUnUn] = conditionalGaussian(mu,Sig,f,varargin)
%% conditionalGaussian
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'mu')
addRequired(Parser,'Sig')
addRequired(Parser,'f')
addParameter(Parser,'nearSPD',true)
addParameter(Parser,'x_indices',NaN)
addParameter(Parser,'SigUnUn',NaN)
addParameter(Parser,'SigObsObs',NaN)
addParameter(Parser,'SigUnObs',NaN)
addParameter(Parser,'SigObsUn',NaN)

parse(Parser,mu,Sig,f,varargin{:})

mu = Parser.Results.mu;
Sig = Parser.Results.Sig;
f = Parser.Results.f;
nearSPD = Parser.Results.nearSPD;
x_indices = Parser.Results.x_indices;
SigUnUn = Parser.Results.SigUnUn;
SigObsObs = Parser.Results.SigObsObs;
SigUnObs = Parser.Results.SigUnObs;
SigObsUn = Parser.Results.SigObsUn;

if isnan(x_indices)
    x_indices = false(size(mu));
    x_indices(1:length(f)) = true;     % Assume observations are of first q elements of length N vector mu
end

N = length(mu);
q = sum(x_indices);

%% Reshape mean and covariance for partitioning
if any(isnan(SigUnUn(:))) || any(isnan(SigObsObs(:))) || any(isnan(SigUnObs(:)))
    indsOrig = 1:N;
    indsNew = [indsOrig(~x_indices) indsOrig(x_indices)];
    % [Is,Js] = meshgrid(indsOrig);
    [Inew,Jnew] = meshgrid(indsNew);
    SigNew = nan(size(Sig));
    for i = 1:N
        for j = 1:N
            SigNew(i,j) = Sig(Inew(i,j),Jnew(i,j));     % Now covariance matrix is of form SigNew = [SigUnUn, SigUnObs; SigObsUn, SigObsObs];
        end
    end
    SigUnUn = SigNew(1:(N-q),1:(N-q));
    SigUnObs = SigNew(1:(N-q),(N-q+1):end);
    SigObsObs = SigNew((N-q+1):end,(N-q+1):end);
    SigObsUn = SigNew((N-q+1):end,1:(N-q));
end

%% Partition mean and covariance
muObs = mu(x_indices);
muUn = mu(~x_indices);


mu_ = muUn + SigUnObs/SigObsObs*(f-muObs);
Sig_ = SigUnUn - SigUnObs/SigObsObs*SigObsUn;
if nearSPD && ~any(isnan(Sig_(:)))
    Sig_ = nearestSPD(Sig_);
end

