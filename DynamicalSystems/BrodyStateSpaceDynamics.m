function [du, dv, ncs] = BrodyStateSpaceDynamics(us,vs,varargin)
%% BrodyNetwork
%
%   [u, v, du, dv] = BrodyNetwork(uinit,vinit)
%
%   Numerically simulates a two neuron model in the vein of Machens et. al.
%   2005 in the "decision" regime.
%
%
%%

%% Defaults
f_default = @(x)(ioMachensBrody(x));
xgrd_default = linspace(-0.0049,0.0075,1234);

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'us')
addRequired(Parser,'vs')
addParameter(Parser,'t',0:1000)
addParameter(Parser,'tau',80)
addParameter(Parser,'f',f_default)
addParameter(Parser,'lambda',2.5625)
addParameter(Parser,'gE',0.0015)
addParameter(Parser,'wI',0.0012)
addParameter(Parser,'xgrd',xgrd_default)
addParameter(Parser,'plotflg',false)

parse(Parser,us,vs,varargin{:})

us = Parser.Results.us;
vs = Parser.Results.vs;
t = Parser.Results.t;
tau = Parser.Results.tau;
f = Parser.Results.f;
lambda = Parser.Results.lambda;
gE = Parser.Results.gE;
wI = Parser.Results.wI;
xgrd = Parser.Results.xgrd;
plotflg = Parser.Results.plotflg;

%% Run simulation
du = -us + f( -wI*vs + lambda*gE );
dv = -vs + f( -wI*us + lambda*(gE + 0.0001) );

% du = nan(length(us),1);
% dv = nan(length(us),1);
% 
% for ui=1:length(us);
%     du(ui) = -us(ui) + f( -wI*vs(ui) + lambda*gE );
%     dv(ui) = -vs(ui) + f( -wI*us(ui) + lambda*gE );
% end

%% Calculate nullclines
ncs.u = 1000*wI*ioMachensBrody( -xgrd + lambda*gE );
ncs.v = 1000*wI*ioMachensBrody( -xgrd + lambda*(gE + 0.0001) );
ncs.x = 1000*xgrd;

%% Plotting
if plotflg
    
%     % Plot trajectories in state space
%     figure('Name','State space trajectory')
%     plot(u,v,'k.-')
%     hold on
%     plot(1000*xgrd,ncs.v,...
%         'Color',[0.5 0.5 0.5])
%     plot(wI*ncs.u,1000*xgrd,...
%         'Color',[0.5 0.5 0.5])
%     axis square
%     axis([0 max(1000*xgrd) 0 max(1000*xgrd)])
%     xlabel('u')
%     ylabel('v')
%     
%     % Plot trajectories over time
%     figure('Name','u and v over time')
%     subplot(1,2,1)
%     plot(t,u,'k')
%     xlabel('t');
%     ylabel('u');
%     subplot(1,2,2)
%     plot(t,v,'k')
%     xlabel('t')
%     ylabel('v')
end

%% Functions
function y = ioMachensBrody ( xin )

initpar;
xgrd = par.x;
fx = par.fx;

N     = length( xgrd );
xmax  = max( xgrd );
fxmax = max( fx );

y  = zeros( size(xin) );
u1 = find( xmax <=xin );
u2 = find( 0<xin & xin<xmax );
u3 = find( xin<=0 );

if ~isempty(u1); y(u1) = fxmax; end
if ~isempty(u2);
  ind = 1+round( N * (xin(u2)-xgrd(1))/ ( xgrd(end)-xgrd(1) ));
  ind( ind>N ) = N;
  y(u2) = fx( ind );
end
if ~isempty(u3); y(u3) = 0; end;