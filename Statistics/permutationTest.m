function p = permutationTest(x,y,varargin)
%% permutationTest
%
%   p = permutationTest(x,y)
%
%%

%% Defaults

%% Parse inputs

Parser = inputParser;

addRequired(Parser,'x')                 % Input data with "x" label
addRequired(Parser,'y')                 % Input data with "y" label
addParameter(Parser,'tail','right')     % Defaults to testing if x is greater than y
addParameter(Parser,'Permutations',NaN) % Number of permutations to perform to generate test distribution; defaults to all permutations

parse(Parser,x,y,varargin{:})

x = Parser.Results.x;
y = Parser.Results.y;
tail = Parser.Results.tail;
Permutations = Parser.Results.Permutations;

%% Find the difference in the means of the data
xbar = mean(x(:));
ybar = mean(y(:));

d = xbar-ybar;

%% Permute the data and find the distribution of d under the null hypothesis
nx = numel(x);
ny = numel(y);

all = [x(:); y(:)];
if isnan(Permutations)
    ALL = perms(all);
else
    for i = 1:Permutations
        inds = randperm(numel(all));
        ALL(i,:) = all(inds);
    end
end
XBars = mean(ALL(:,1:nx),2);
YBars = mean(ALL(:,nx+1:ny),2);
D = XBars-YBars;

switch tail
    case {'left'}
        p = sum(D < d)/numel(D);
        
    case {'right'}
        p = sum(D > d)/numel(D);
        
    otherwise
        error(['tail option ' tail ' not recognized!'])
end