function Y = Isomap(X,outputDim,varargin)
%% Isomap
%
%   Y = Isomap(X,outputDim)
%   Computes the isomap (Tenenbaum et. al. 2000) of the data (X) on a lower
%   dimensional manifold (Y) of dimensionality outputDim based on the
%   Euclidean distance of the points from one another.
%
%   Y = Isomap(...,'neighborRule',neighborRule)
%   Computes the isomap for neighbors specified by the rule neighborRule
%       neighborRule.type = 'K-closest' - K closest neighbors (default)
%                           with K specified by neightborRule.K (default =
%                           10).
%       neightborRule.type = 'WithinEpsilon' - Points that are within
%                            neighborRule.Ellispon uits (Euclidean
%                            distance).
%
%   Y = Isomap(...,'distances',D)
%   Allows the user to provide a set of "distances" made according to
%   his/her own distance metric.
%
% Notes: based on the Isomap algorithm by Tenebaum et. al. (2000) "A global
%        geometric framework for nonlinear dimensionality reduction."
% 
%   Date:       author:    Comment:
%   2015.4.8    swe        Initial commit
%
%%

% Defaults
neighborRule_default.type = 'K-closest';
neighborRule_default.K = 10;

% Parse Inputs
Parser = inputParser;

addRequired(Parser,'X');
addRequired(Parser,'outputDim');
addParameter(Parser,'neighborRule',neighborRule_default);
addParameter(Parser,'distances',NaN);
addParameter(Parser,'TestY',NaN);

parse(Parser,X,outputDim,varargin{:})

X = Parser.Results.X;
outputDim = Parser.Results.outputDim;
neighborRule = Parser.Results.neighborRule;
distances = Parser.Results.distances;
TestY = Parser.Results.TestY;


if isnan(distances)
    % Find Euclidean distances between points if distance measurements are not
    % provided
    D = permute(sqrt(sum((repmat(X,[1 1 size(X,1)]) - permute(repmat(X,[1 1 size(X,1)]),[3 2 1])).^2,2)),[1 3 2]);
end
INF =  1000*max(max(D))*size(D,1);

% Find the neighbors on the lower dimensional manifold
switch neighborRule.type
    case 'K-closest'
%         [~, indx] = sort(D,2);
%         DG = Inf(size(D));
%         for i = 1:size(D,1)
%             DG(i,indx(i,1:neighborRule.K)) = D(i,indx(i,1:neighborRule.K));
%         end
        
        for i = 1:size(D,2)
            [a, indx] = sort(D(:,i));
            for j = 1:size(D,1)
                if D(j,i) <= a(neighborRule.K)
                    DG(j,i) = D(j,i);
                else
                    DG(j,i) = INF;
                end
            end
        end
        
    case 'WithinEpsilon'
        error('epsilon-Isomap not yet implemented!')
        
end

% Compute the minimum path length along the graph from each point to all other points
% for i = 1:size(D,1)
%     for j = 1:size(D,2)
%         for k = 1:size(D,1)
%             %disp([k DG(i,j) DG(i,k) DG(k,j) DG(i,k)+DG(k,j)])
%             DG(i,j) = min([DG(i,j) DG(i,k)+DG(k,j)]);
%         end
%     end
% end
% 
% [rowM, colM] = meshgrid(1:size(D,1));
% for i = 1:size(D,1)
%     r = rowM == i;
%     c = colM == i;
%     Dtemp = DG*r + DG*c;
%     Dtemp(Dtemp == 0) = Inf;
%     DG = min(cat(3,DG,Dtemp),[],3);
% end

for k=1:size(D,1)
     DG = min(DG,repmat(DG(:,k),[1 size(D,1)])+repmat(DG(k,:),[size(D,1) 1])); 
end

% Construct d-dimensional embedding
taus = tau(DG);
[v, d] = eig(taus);
% d = real(diag(d));
% v = real(v);
% Y = zeros(size(X,1),outputDim);
% 
% for i = 1:outputDim
%     Y(:,i) = real(sqrt(d(end-i+1))*v(:,end-i+1));
% end

h = real(diag(d)); 
[~, sorth] = sort(h);
sorth = sorth(end:-1:1); 
d = real(diag(d(sorth,sorth))); 
v = real(v(:,sorth)); 

for i = 1:outputDim
    Y(:,i) = real(sqrt(d(i))*v(:,i));
end

end
%% Functions
function out = tau(in)
S = in.^2;
H = eye(size(in)) - 1/size(in,1);
out = -H*S*H/2;
end