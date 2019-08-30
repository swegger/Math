function X = CostToGo(t,Q,R,A,B,varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'t')
addRequired(Parser,'Q')
addRequired(Parser,'R')
addRequired(Parser,'A')
addRequired(Parser,'B')

parse(Parser,t,Q,R,A,B,varargin{:})

t = Parser.Results.t;
Q = Parser.Results.Q;
R = Parser.Results.R;
A = Parser.Results.A;
B = Parser.Results.B;

%% Discrete-time algebreic Riccati equation
X = nan(size(Q,1),size(Q,2),length(t));
X(:,:,end) = Q;
for i = 1:length(t)-1
    ti = length(t)-i;
    X(:,:,ti) = Q + A'*X(:,:,ti+1)*A -...
        A'*X(:,:,ti+1)*B*inv(B'*X(:,:,ti+1)*B + R)*B'*X(:,:,ti+1)*A;
end
