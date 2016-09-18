function [u, w, du, dw] = ExpExpLinearModel(t,u0,w0,A,B,I)
%% SigExpLinearModel
%
%   [u, w, du, dw] = SigExpLinearModel(t,u0,w0,A,F,B,I)
%
%   Simulates a nonlinear-linear dynamic model of the form:
%       du = dt * ( A(1,1)*( u/(1+u^2) ) + A(1,2)*( exp(w) ) + I(1) + eta(1)
%       dw = dt * ( A(2,1)*( exp(w) ) + A(2,1) * ( u/(1+u^2) ) + I(2) + eta(2)
%
%%

if size(I,2) == 1
    I = repmat(I,1,length(t));
end

dt = t(2)-t(1);
du = nan(size(dt));
dw = nan(size(dt));
for ti = 1:length(t)
    if ti == 1
        eta = B*randn(2,1);
        du(ti) = ( A(1,1)*( exp(u0) ) +...
            A(1,2)*( exp(w0) ) + I(1,ti) + eta(1) )*dt;
        dw(ti) = ( A(2,1)*( exp(u0) ) +...
            A(2,2)*( exp(w0) ) + I(2,ti) + eta(2) )*dt;
%         du(ti) = ( A(1,1)*( u0 ) +...
%             A(1,2)*( exp(w0) ) + I(1,ti) + eta(1) )*dt;
%         dw(ti) = ( A(2,1)*( u0 ) +...
%             A(2,2)*( exp(w0) ) + I(2,ti) + eta(2) )*dt;
        u(ti) = u0;
        w(ti) = w0;
    else
        eta = B*randn(2,1);
        du(ti) = ( A(1,1)*( exp(u(ti-1)) ) +...
            A(1,2)*( exp(w(ti-1)) ) + I(1,ti) + eta(1) )*dt;
        dw(ti) = ( A(2,1)*( exp(u(ti-1)) ) +...
            A(2,2)*( exp(w(ti-1)) ) + I(2,ti) + eta(2) )*dt;
%         du(ti) = ( A(1,1)*( u(ti-1) ) +...
%             A(1,2)*( exp(w(ti-1)) ) + I(1,ti) + eta(1) )*dt;
%         dw(ti) = ( A(2,1)*( u(ti-1) ) +...
%             A(2,2)*( exp(w(ti-1)) ) + I(2,ti) + eta(2) )*dt;
        u(ti) = u(ti-1) + du(ti);
        w(ti) = w(ti-1) + dw(ti);
    end
end