% 
% FHN matlab code
% J Rinzel 2006
% FHN eqns:
% 
% dv/dt = -v*(v-a)*(v-1)-w+I0
% dw/dt = eps*(v-d*w)
% 
% This code also plots nullclines of the FHN eqns.

% +++++++++++++++++++

function FHN
% dimensionless model of neuronal excitability
% v=Y(1); w=Y(2);
% a=0.1;d=1;eps=0.01;I0=0; v0=-0.25;w0=0;%no input, no IC
a=0.1;d=1;eps=0.01;I0=.1; v0=0.;w0=0;%Some input, no IC
% a=0.1; d=1;eps=0.01;I0=0; v0=-0.1;w0=-0.1; %No input, some IC
% a=0.1; d=1;eps=0.01;I0=0; v0=0.25;w0=0; %No input, some IC
Y0=[v0,w0];
t=0:0.1:400;
options=odeset('RelTol',1.e-5);
[T, Y]=ode45(@dydt_FHN,t,Y0,options,a,eps,d,I0);
figure('Position',[177   380   290   214])
plot(T,Y(:,1),T,Y(:,2)); % time courses of V and w
box off
legend('v(t)','w(t)');
xlabel('Time'); ylabel('v, w');
set(gca,'YLim',[-0.4 1.2])
vpts=(-1.5:.05:1.5);
drawnow();
set(gcf,'Position',[177   380   290   214])
drawnow();
figure('Position',[177   380   290   214])
hold on;
plot(Y(:,1),Y(:,2)); % V-w phase plane 
line(Y(1,1),Y(1,2),'Marker','o')
line(Y(end,1),Y(end,2),'Marker','x')
%determine and plot the v,w-nullclines
%options=optimset; % sets options in fzero to default values
%for k=1:61
%vnullpts(k)=fzero(@vrhs_FHN,[-10 10],options,vpts(k),a,I0);  
%end
vnullpts=-vpts.*(vpts-a).*(vpts-1)+I0;
wnullpts=vpts/d;
plot(vpts,vnullpts,'Color',[0.5 0.5 0.5])
plot(vpts,wnullpts,'Color',[0.5 0.5 0.5]);
xlabel('v'); ylabel('w');
axis([-0.5 1.25 -0.1 0.45]);
drawnow();
set(gcf,'Position',[177   380   290   214])
drawnow();

% +++++++++++++++++++

function dY=dydt_FHN(t,Y,a,eps,d,I0)
v=Y(1);
w=Y(2);
dY=zeros(2,1);
% dY(1)=-v*(v-a)*(v-1)-w+I0*1/(1+exp(20-t)/.2); %dunno what the exp does
dY(1)=-v*(v-a)*(v-1)-w+I0;
dY(2)=eps*(v-d*w);

% ++++++++++++++++++++

function val=vrhs_FHN(w,v,a,I0)
	val=-v*(v-a)*(v-1)-w+I0;