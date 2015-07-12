function p = GaussMix(w,mu,sig,x,varargin)
%% GaussMix
%
%   Gaussian mixture model of a probability distribution.
%
%   p = GaussMix(w,mu,sig,x)
%       w - weighting associated with each Gaussian; 1xK vector, K = # of
%           Gaussians
%       mu - mean of each Gaussian; KxD matrix, D = # of dimensions
%       sig - Covariance of each Gaussian; DxDxK matrix
%       x - points to evaluate probability; NxD matrix, N = # of points
%
%   by Seth Egger
%       Update       Author       Note
%       2014/05/07   swe          Commit with comments
%
%%

% Parse inputs
Parser = inputParser;
addRequired(Parser,'w');
addRequired(Parser,'mu');
addRequired(Parser,'sig');
addRequired(Parser,'x');
addParameter(Parser,'Batch',0);

parse(Parser,w,mu,sig,x,varargin{:})

w = Parser.Results.w;
mu = Parser.Results.mu;
sig = Parser.Results.sig;
x = Parser.Results.x;
Batch = Parser.Results.Batch;


% Force w to be a row vector
w = w(:)';

% Number of Gaussians
K = length(w);

% Number of Dimensions
D = size(x,2);

% Number of evaluations
if Batch == 0
    N = size(x,1);
    
    % Set up variables
    MU = reshape(permute(repmat(mu,[1 1 N]),[3 1 2]),[N*K,D,1]);
    SIG = reshape(permute(repmat(sig,[N,1,1]),[2 1 3]),[D D N*K]);
    X = repmat(x,K,1);
    
    P = mvnpdf(X,MU,SIG);
    
    p = (w*reshape(P,[N,K])')';
else
    for i = 1:ceil(size(x,1)/Batch)
        if length(x) <= (i-1)*Batch+Batch
            xtemp = x((i-1)*Batch+1:end,:);
            N = size(xtemp,1);
            
            % Set up variables
            MU = reshape(permute(repmat(mu,[1 1 N]),[3 1 2]),[N*K,D,1]);
            SIG = reshape(permute(repmat(sig,[N,1,1]),[2 1 3]),[D D N*K]);
            X = repmat(xtemp,K,1);
            
            P = mvnpdf(X,MU,SIG);
            
            p((i-1)*Batch+1:size(x,1),:) = (w*reshape(P,[N,K])')';
        else
            xtemp = x((i-1)*Batch+1:(i-1)*Batch+Batch,:);
            N = size(xtemp,1);
            
            % Set up variables
            MU = reshape(permute(repmat(mu,[1 1 N]),[3 1 2]),[N*K,D,1]);
            SIG = reshape(permute(repmat(sig,[N,1,1]),[2 1 3]),[D D N*K]);
            X = repmat(xtemp,K,1);
            
            P = mvnpdf(X,MU,SIG);
            
            p((i-1)*Batch+1:(i-1)*Batch+Batch,:) = (w*reshape(P,[N,K])')';
        end
    end
end