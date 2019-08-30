function u = optimizedControl(t,y0,Q,R,A,B)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'t')
addRequired(Parser,'y0')
addRequired(Parser,'Q')
addRequired(Parser,'R')
addRequired(Parser,'A')
addRequired(Parser,'B')

parse(Parser,t,y0,Q,R,A,B)

t = Parser.Results.t;
y0 = Parser.Results.y0;
Q = Parser.Results.Q;
R = Parser.Results.R;
A = Parser.Results.A;
B = Parser.Results.B;

%% Backwards induction of control solution
u = nan(size(B,2),length(t));

X = CostToGo(t,Q,R,A,B);

y(:,1) = y0;
u(:,1) = -inv(R + B'*X(:,:,1)*B)*(B'*X(:,:,1)*A)*y(:,1);
for ti = 2:length(t)-1
    y(:,ti) = A*y(:,ti-1) + B*u(ti-1);
    u(:,ti) = -inv(R + B'*X(:,:,ti+1)*B)*(B'*X(:,:,ti+1)*A)*y(:,ti);
end