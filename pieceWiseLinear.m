function tout = pieceWiseLinear(tin,p)
%% pieceWise
%
%
%%
m = p(1);
b = p(2);
t0 = p(3);
tout = zeros(size(tin));

tout(tin >= p(3)) = m*tin( tin >= p(3) ) + b;

tout(tin < p(3)) = m*t0 + b;

