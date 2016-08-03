function [x,rho,log] = rqmin1(A,M,x,tol,C)
%RQMIN1 [x,rho] = rqmin1(A,M,x0,tol,C)
% cg-Rayleigh-Quotienten-Minimization for the computation
% of the smallest eigenvalue of A*x = lambda*M*x,
% A and M are symmetric, M spd. x0 initial vector
% C?*C preconditioner
% tol: convergence criterium:
% ||2*(C?*C)\(A*x - lam*M*x)|| < tol
% PA 16.6.2000
u = M*x;
q = sqrt(x'*u);
x = x/q; u = u/q;
v = A*x;
rho = x'*v;
k = 0; g = x; gnorm = 1; log=[]; % Initialisierungen
while gnorm > tol,
k = k + 1;
galt = g;
if exist('C'),
g = 2*(C\(C'\(v - rho*u))); % vorkonditionierter Gradient
else
g = 2*(v - rho*u); % Gradient
end
if k == 1,
p = -g;
else
p = -g + (g'*M*g)/(galt'*M*galt)*p;
end
[qq,ll] = eig([x p]'*[v A*p],[x p]'*[u M*p]);
[rho,ii] = min(diag(ll));
delta = qq(2,ii)/qq(1,ii);
x = x + delta*p;
u = M*x;
q = sqrt(x'*u);
x = x/q; u = u/q;
v = A*x;
gnorm = norm(g);
if nargout>2, log = [log; [k,rho,gnorm]]; end
end