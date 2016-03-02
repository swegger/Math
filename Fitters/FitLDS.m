 function [w, mu, sig, logl, gamma] = FitLDS(X,mu0init,V0init,Ainit,Ginit,Cinit,Sinit,tolerance,varargin)
%% FitMixtureModelEM
%
%   [w, mu, sig, logl, gamma] = FitLDS(x,Ainit,Ginit,Cinit,Sinit,tolerance)
%   
%   Fits a Gaussina linear dynamic system of the form:
%       z(k+1) = Az(k) + w(k)
%       x(k) = Cz(k) + v(k)
%       w(k) ~ N(w|0,G)
%       v(k) ~ N(v|0,S)
%
%%

%% Defaults
noiseModel_default.type = 'none';

%% Parse inputs
Parser = inputParser;
addRequired(Parser,'X');
addRequired(Parser,'Ainit');
addRequired(Parser,'Ginit');
addRequired(Parser,'Cinit');
addRequired(Parser,'Sinit');
addRequired(Parser,'tolerance');
addParameter(Parser,'Batch',0);
addParameter(Parser,'percentFit',100);
addParameter(Parser,'verbose','logLikelihood');
addParameter(Parser,'noiseModel',noiseModel_default);
addParameter(Parser,'method','EM');

parse(Parser,X,Ainit,Ginit,Cinit,Sinit,tolerance,varargin{:})

X = Parser.Results.X;
Ainit = Parser.Results.Ainit;
Ginit = Parser.Results.Ginit;
Cinit = Parser.Results.Cinit;
Sinit = Parser.Results.Sinit;
tolerance = Parser.Results.tolerance;
Batch = Parser.Results.Batch;
percentFit = Parser.Results.percentFit;
verbose = Parser.Results.verbose;
noiseModel = Parser.Results.noiseModel;
method = Parser.Results.method;

%% Fit the model to the data
switch method
    case 'EM'
        % Find Z based on initial conditions
        Zhat = nan(size(X,1),length(mu0init));
        Zhat(1,:) = mu0init;
        for i = 2:size(X,1)
            Zhat(i,:) = Zhat(i-1) + Ainit*Zhat(i-1,:)';
        end
        dif = 10^6;
        while dif >= tolerance | difLL >= lltolerance
            % Find the new estimate of Z, given the parameters
            Zhat(1,:) = mu0;
            for i = 2:size(X,1)
                Zhat(i,:) = Zhat(i-1) + A*Zhat(i-1,:)';
            end
            
            % Find the new initial mean and covariance
            mu0 = ;
            V0 = ;
            
            % Find the new transition parameters
            A = 
        
    otherwise
        error(['Fitting method ' method ' not recognized!'])
end