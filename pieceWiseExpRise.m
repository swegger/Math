function R = pieceWiseExpRise(tin,p)
%% pieceWise
%
%
%%
R0 = p(1);
Rf = p(2);
tau = p(3);
t0 = p(4);
R = zeros(size(tin));

R(tin >= t0) = R0 + (Rf-R0)*(1-exp(-(tin(tin>=t0)-t0)/tau));

%R(tin < t0,:) = R0;
R(tin < t0) = R0;

