function [P, feval, exitflg, output] = Fit_Lisberger1999(t,v,varargin)
%% Fit_Lisberger1999
%
%   P = Fit_Lisberger1999(t,v)
%
%   Fits the parameters of the function of eye velocity over time from the
%   onset of target motion:
%
%       v(1) = P(4)*(t-P(1))^P(2) / ( (t-P(1))^P(2) + P(3)^P(2) )
%       v(2) = P(5)*(t-P(1))^P(2) / ( (t-P(1))^P(2) + P(3)^P(2) )
%
%%

%% Defaults
OPTIONS = optimset('Display','iter');

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'t')
addRequired(Parser,'v')
addParameter(Parser,'P0',[50, 1, 10, 5, 15])
addParameter(Parser,'options',OPTIONS)
addParameter(Parser,'LB',[0,0,0,0,0])
addParameter(Parser,'UB',[250,Inf,Inf,Inf,Inf])

parse(Parser,t,v,varargin{:})

t = Parser.Results.t;
v = Parser.Results.v;
P0 = Parser.Results.P0;
LB = Parser.Results.LB;
UB = Parser.Results.UB;
options = Parser.Results.options;

%% Fit the data
f = @(p)(minimizant(p,t,v));
[P, feval, exitflg, output] = fmincon(f, P0, [], [], [], [],LB,UB,[],options);

%% Functions

%% Minimizant
function error = minimizant(p,t,v)
fnc = Lisberger1999(t,p);
error = nansum((v(:) - fnc(:)).^2);

%% Lisberger1999
function out = Lisberger1999(t,p)
out(:,1) = p(4)*(t(:,1)-p(1)).^p(2) ./ ( (t(:,1)-p(1)).^p(2) + p(3)^p(2) );
out(:,2) = p(5)*(t(:,2)-p(1)).^p(2) ./ ( (t(:,2)-p(1)).^p(2) + p(3)^p(2) );
out(t(:,1) < p(1),1) = 0;
out(t(:,2) < p(1),2) = 0;
% out(isnan(out)) = 0;