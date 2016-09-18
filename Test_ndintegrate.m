function Test_ndintegrate
%% Test_ndintegrate
%
%   Tests ndintegrate with various user defined funcitons:
%
%       f(x,y) = y(1)*x(1) + y(2)*x(2) + y(3);
%           N = 2; (Number of dimensions to set below)
%       f2(x,y) = y(1)*x(1)^2 + y(2)*x(2)^2 + y(3);
%           N = 2;
%       f3(x,y) = y(1)*x(1) + y(2)*x(2) + y(3)*x(3);
%           N = 3;
%
%   Designed to have a breakpoint at "stopline" so that the user can review
%   the results of the integration.
%
%%

% Set up problem
N = 2;
MinMax = repmat([0 1],N,1);
funcitonHandle = @(x,y)(f2(x,y));
options.dx = 0.01;          % dx = 1% of the distance from min to max for quad and trapz
options.N = 10000;          % Number of samples to use in Monte Carlo integration
Y = [1 1 1; 1 2 3; 1 0 0; 0 0 1];

tic
[integraldx(:,1), errest{1}] = ndintegrate(funcitonHandle,MinMax,'method','quad','options',options,'ExtraVariables',Y);
toc
tic
[integraldx(:,2), errest{2}] = ndintegrate(funcitonHandle,MinMax,'method','trapz','options',options,'ExtraVariables',Y);
toc
tic
[integraldx(:,3), errest{3}] = ndintegrate(funcitonHandle,MinMax,'method','MonteCarlo','options',options,'ExtraVariables',Y);
toc

stopline

%% Functions
% Linear combination of x using parameters in y: f(x,y) = y(1)*x(1) + y(2)*x(2) + y(3)
function out = f(x,y)
% Arrange variables so that each row corresponds to a different combination
% of x's, each column a different combination of y's and the third
% dimension are the dimensions of x
X = repmat(permute([x ones(size(x,1),1)],[1 3 2]),1,size(y,1));
Y = repmat(permute(y,[3 1 2]),size(X,1),1);

% Sum across x(1) and x(2)
out = sum(X.*Y,3);


% Linear combination of x^2 using parameters in y: f2(x,y) = y(1)*x(1)^2 + y(2)*x(2)^2 + y(3)
function out = f2(x,y)
% Arrange variables so that each row corresponds to a different combination
% of x's, each column a different combination of y's and the third
% dimension are the dimensions of x
X = repmat(permute([x ones(size(x,1),1)],[1 3 2]),1,size(y,1));
Y = repmat(permute(y,[3 1 2]),size(X,1),1);

% Sum across x(1) and x(2)
out = sum(X.^2.*Y,3);

% Linear combination of x using parameters in y: f3(x,y) = y(1)*x(1) + y(2)*x(2) + y(3)*x(3)
function out = f3(x,y)
% Arrange variables so that each row corresponds to a different combination
% of x's, each column a different combination of y's and the third
% dimension are the dimensions of x
X = repmat(permute(x,[1 3 2]),1,size(y,1));
Y = repmat(permute(y,[3 1 2]),size(X,1),1);

% Sum across x(1) and x(2)
out = sum(X.*Y,3);