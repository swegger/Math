% LLE ALGORITHM (using K nearest neighbors)
%
% [Y] = lle(X,K,dmax)
%
% X = data as D x N matrix (D = dimensionality, N = #points)
% K = number of neighbors
% dmax = max embedding dimensionality
% Y = embedding as dmax x N matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y] = myLLE(X,K,d,varargin)

%% Defaults
LabelOpts_default.Labels = NaN;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'X');
addRequired(Parser,'K');
addRequired(Parser,'d');
addParameter(Parser,'LabelOpts',LabelOpts_default)

parse(Parser,X,K,d,varargin{:})

X = Parser.Results.X;
K = Parser.Results.K;
d = Parser.Results.d;
LabelOpts = Parser.Results.LabelOpts;

%% Handle data labels
% Validate LabelOpts
if isnan(LabelOpts.Labels)
    LabelOpts.Labels = ones(1,size(X,2));
    LabelOpts.minPercent = 1;
elseif ~isfield(LabelOpts,'minPercent')
    LabelOpts.minPercent = 1;
end
L = unique(LabelOpts.Labels);
nControl = ceil(K*LabelOpts.minPercent);
nper = floor(nControl/length(L));

%% Perform LLE

[D,N] = size(X);
fprintf(1,'LLE running on %d points in %d dimensions\n',N,D);


% STEP1: COMPUTE PAIRWISE DISTANCES & FIND NEIGHBORS 
fprintf(1,'-->Finding %d nearest neighbours.\n',K);

X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;

[~,index] = sort(distance);
neighborhood = nan(K,N);
if length(L) > 1
    for ni = 1:N
        ntemp = nan(K,1);
        for li = 1:length(L)
            indtemp = index(LabelOpts.Labels' == L(li) & index(:,ni) ~= index(1),ni);
            ntemp(nper*(li-1)+1:nper*(li-1)+nper) = indtemp(1:nper);
        end
        remainingSpots = sum(isnan(ntemp));
        startind = 2;
        while remainingSpots > 0
            newSpots = setdiff(index(startind:startind+remainingSpots-1,ni),...
                ntemp(~isnan(ntemp)));
            ntemp(end-remainingSpots+1:end-remainingSpots+length(newSpots)) = newSpots;
            remainingSpots = sum(isnan(ntemp));
            startind = startind+remainingSpots;
        end
        
        neighborhood(:,ni) = sort(ntemp(1:K));
    end
else
    neighborhood = index(2:(1+K),:);
end
    


% STEP2: SOLVE FOR RECONSTRUCTION WEIGHTS
fprintf(1,'-->Solving for reconstruction weights.\n');

if(K>D) 
  fprintf(1,'   [note: K>D; regularization will be used]\n'); 
  tol=1e-3; % regularlizer in case constrained fits are ill conditioned
else
  tol=0;
end

W = zeros(K,N);
for ii=1:N
   z = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K); % shift ith pt to origin
   C = z'*z;                                        % local covariance
   C = C + eye(K,K)*tol*trace(C);                   % regularlization (K>D)
   W(:,ii) = C\ones(K,1);                           % solve Cw=1
   W(:,ii) = W(:,ii)/sum(W(:,ii));                  % enforce sum(w)=1
end;


% STEP 3: COMPUTE EMBEDDING FROM EIGENVECTS OF COST MATRIX M=(I-W)'(I-W)
fprintf(1,'-->Computing embedding.\n');

% M=eye(N,N); % use a sparse matrix with storage for 4KN nonzero elements
M = sparse(1:N,1:N,ones(1,N),N,N,4*K*N); 
for ii=1:N
   w = W(:,ii);
   jj = neighborhood(:,ii);
   M(ii,jj) = M(ii,jj) - w';
   M(jj,ii) = M(jj,ii) - w;
   M(jj,jj) = M(jj,jj) + w*w';
end;

% CALCULATION OF EMBEDDING
options.disp = 0; options.isreal = 1; options.issym = 1; 
[Y,eigenvals] = eigs(M,d+1,0,options);
Y = Y(:,2:d+1)'*sqrt(N); % bottom evect is [1,1,1,1...] with eval 0


fprintf(1,'Done.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% other possible regularizers for K>D
%   C = C + tol*diag(diag(C));                       % regularlization
%   C = C + eye(K,K)*tol*trace(C)*K;                 % regularlization
