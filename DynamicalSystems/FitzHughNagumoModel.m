function [u, w, du, dw] = FitzHughNagumoModel(t,u0,w0,epsilon,b,I,a)
%% FitzHughNagumoModel
%
%   [u, w] = FizHugNagumoModel(t,u0,w0,epsilon,b,I)
%
%   Implements the FitzHugh-Nagmomo 2D dynamical system model.
%
%%

if length(epsilon) == 1
    epsilon = epsilon*ones(size(t));
end
if size(b,1) == 1
    b = repmat(b,length(t),1);
end
if length(I) == 1
    I = I*ones(size(t));
end

dt = t(2)-t(1);
for ti = 1:length(t)
    if ti == 1
        du(ti) = (u0 - u0.^3/3 - w0 + I(ti) + a*randn)*dt;
        dw(ti) = (epsilon(ti)*(b(ti,1) + b(ti,2)*u0 - w0) + a*randn)*dt;
        u(ti) = u0;
        w(ti) = w0;
    else
        du(ti) = (u(ti-1) - u(ti-1).^3/3 - w(ti-1) + I(ti) + a*randn)*dt;
        dw(ti) = (epsilon(ti)*(b(ti,1) + b(ti,2)*u(ti-1) - w(ti-1)) + a*randn)*dt;
        u(ti) = u(ti-1) + du(ti);
        w(ti) = w(ti-1) + dw(ti);
    end
end
        