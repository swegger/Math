function [du, dw, u_null, w_null, us] = ExpExpLinearModelVectorField(u,w,A,I)
%% FitzHughNagumoModelVectorField
%
%
%%

% Vector field
du = A(1,1)*( exp(u) ) + A(1,2)*( exp(w) );% + I(1);
dw = A(2,1)*( exp(u) ) + A(2,2)*( exp(w) );% + I(2);

% du = A(1,1)*( u ) + A(1,2)*( exp(w) );% + I(1);
% dw = A(2,1)*( u ) + A(2,2)*( exp(w) );% + I(2);

% Nullclines
if isdiag(A)
    us = unique(u);
    u_null = sqrt( I(1) / (I(1) -1) );
    w_null = log( I(2) );
else
    u_null = NaN;
    w_null = NaN;
end
