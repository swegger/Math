function [x, dx] = WKWmodel(x0,W,K,Gamma,epsilon,varargin)
%% WKWmodel
%
%   [x,dx] = WKWmodel(x0,W,K,Gamma)
%
%   Implements WKWmodel variation of two neuron model
%
%%

%% Defaults
f_default = @(x)(ioMachensBrody(x));
xgrd_default = linspace(-0.0049,0.0075,1234);

%% Input parser
Parser = inputParser;

addRequired(Parser,'x0')            % Initial conditions
addRequired(Parser,'W')
addRequired(Parser,'K')
addRequired(Parser,'Gamma')
addRequired(Parser,'epsilon')
addParameter(Parser,'t',0:1000)
addParameter(Parser,'tau',80)
addParameter(Parser,'f',f_default)
addParameter(Parser,'lambda',2.5625)
addParameter(Parser,'gE',0.0015)
addParameter(Parser,'transitionNoise',0)
addParameter(Parser,'xgrd',xgrd_default)
addParameter(Parser,'inputNoise',0)
addParameter(Parser,'trialN',1);
addParameter(Parser,'plotflg',false)

parse(Parser,x0,W,K,Gamma,epsilon,varargin{:})

x0 = Parser.Results.x0;
W = Parser.Results.W;
K = Parser.Results.K;
Gamma = Parser.Results.Gamma;
t = Parser.Results.t;
tau = Parser.Results.tau;
f = Parser.Results.f;
lambda = Parser.Results.lambda;
gE = Parser.Results.gE;
transitionNoise = Parser.Results.transitionNoise;
xgrd = Parser.Results.xgrd;
inputNoise = Parser.Results.inputNoise;
trialN = Parser.Results.trialN;
plotflg = Parser.Results.plotflg;
epsilon = Parser.Results.epsilon;

t = t(:);
N = size(x0,1);
D = size(K,1);
T = size(t,1);

%% Run simulation
dx = nan([N,trialN,T]);
x = nan([N,trialN,T]);
x(:,:,1) = repmat(x0,[1 trialN 1]);
xlin = nan([N,trialN,T]);
fx = nan([N,trialN,T]);
v = nan([D,trialN,T]);
v(:,:,1) = W'*repmat(x0,[1 trialN 1]);

for ti=2:length(t)
    
    % ODE integration using Euler's method
    I = gE(:)' + inputNoise*randn(1,trialN);
    
    xlin(:,:,ti) = W*(K*W'*x(:,:,ti-1) + lambda*Gamma'*x(:,:,ti-1) +...
        lambda*epsilon);
    for xi = 1:N
        for triali = 1:trialN
            fx(xi,triali,ti) = f( xlin(xi,triali,ti) + lambda*I );
            dx(xi,triali,ti) = -x(xi,triali,ti-1) + ...
                fx(xi,triali,ti) + ...
                transitionNoise*randn;
        end
    end
    x(:,:,ti) = x(:,:,ti-1)+dx(:,:,ti)/tau;
    v(:,:,ti) = W'*x(:,:,ti);
end


%% Calculate nullclines
wI = abs(K(1,2));
ncs.v1 = 1000*wI*ioMachensBrody( -repmat(xgrd,[length(gE),1]) +...
    repmat(lambda*gE(:),[1 length(xgrd)]) );
ncs.v2 = 1000*wI*ioMachensBrody( -repmat(xgrd,[length(gE),1]) +...
    repmat(lambda*(gE(:) + 0.0001),[1 length(xgrd)]) );
ncs.x = 1000*xgrd;

%% Plotting
if plotflg
    
    % Plot trajectories in state space
    figure('Name','State space trajectory')
    plot(1000*wI*squeeze(v(1,1,:)),1000*wI*squeeze(v(2,1,:)),'k.-')
    hold on
    plot(1000*wI*squeeze(v(1,1,1)),1000*wI*squeeze(v(2,1,1)),'ko')
    plot(1000*wI*squeeze(v(1,1,end)),1000*wI*squeeze(v(2,1,end)),'ks')
    plot(1000*xgrd,ncs.v2,...
        'Color',[0.5 0.5 0.5])
    plot(ncs.v1,1000*xgrd,...
        'Color',[0.5 0.5 0.5])
    axis square
    axis([0 max(1000*xgrd) 0 max(1000*xgrd)])
    xlabel('v_1')
    ylabel('v_2')
    
    % Plot trajectories over time
    figure('Name','u and v over time')
    for vi = 1:D
        subplot(D,1,vi)
        plot(t,squeeze(v(vi,1,:)),'k')
        xlabel('t');
        ylabel(['v_' num2str(vi)]);
    end
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
if ~isempty(u2)
  ind = 1+round( N * (xin(u2)-xgrd(1))/ ( xgrd(end)-xgrd(1) ));
  ind( ind>N ) = N;
  y(u2) = fx( ind );
end
if ~isempty(u3); y(u3) = 0; end