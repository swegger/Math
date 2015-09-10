function [du, dw] = FitzHughNagumoModelVectorField(u,w,epsilon,b,I)
%% FitzHughNagumoModelVectorField
%
%
%%

du = u - u.^3/3 - w + I;
dw = epsilon*(b(1) + b(2)*u - w);

