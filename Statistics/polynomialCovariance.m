function C = polynomialCovariance(fparams,gparams,x,varargin)
%% polynomialCovariance
%
%   
%
%%

%% Defaults
nonlinearity_default = @(x)(x);

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'fparams')
addRequired(Parser,'gparams')
addRequired(Parser,'x')
addParameter(Parser,'nonlinearity',nonlinearity_default)

parse(Parser,fparams,gparams,x,varargin{:});

fparams = Parser.Results.fparams;
gparams = Parser.Results.gparams;
x = Parser.Results.x;
nonlinearity = Parser.Results.nonlinearity;

%% Calculate covaraince for each supplied pair
fpoly = polyval(fparams,x(:,1)+x(:,2));
gpoly = polyval(gparams,x(:,1)-x(:,2));

f = nonlinearity(fpoly);
g = nonlinearity(gpoly);

C = f.*g;
