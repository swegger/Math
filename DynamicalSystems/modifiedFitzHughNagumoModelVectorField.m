function [du, dw, u_null, w_null, us] = modifiedFitzHughNagumoModelVectorField(u,w,epsilon,b,I,tau)
%% FitzHughNagumoModelVectorField
%
%
%%

% Vector field
du = (u - u.^3/3 - w + I)/tau;
dw = epsilon*(b(1) + b(2)*u - w)/tau;


% Nullclines
us = unique(u);
u_null = us - us.^3/3 + I;
w_null = b(1) + b(2)*us;
