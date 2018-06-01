function [du, dw, u_null, w_null, us] = FitzHughNagumoModelVectorField2(u,w,epsilon,b,I)
%% FitzHughNagumoModelVectorField
%
%
%%

% Vector field
du = -u.*(u-b(2)).*(u-1) -w + I;
dw = epsilon*(u - b(1)*w);

% Nullclines
us = unique(u);
us = linspace(min(us),max(us),100);

u_null = -us.*(us-b(2)).*(us-1) + I;
w_null = us/b(1);

