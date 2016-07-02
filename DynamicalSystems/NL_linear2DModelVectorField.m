function [du, dw, u_null, w_null, us] = NL_linear2DVectorField(u,w,epsilon,b,I)
%% FitzHughNagumoModelVectorField
%
%
%%

% Vector field
du = u - u.^3/3 - w + I;
dw = epsilon*(b(1) + b(2)*u - w);


% Nullclines
us = unique(u);
u_null = us - us.^3/3 + I;
w_null = b(1) + b(2)*us;
