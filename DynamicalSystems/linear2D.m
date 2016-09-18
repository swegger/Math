function [u, du] = linear2D(t,u0,A,I,B)
%% linear2D
%%

if length(I) == 1
    I = I*ones(size(t));
end

dt = t(2)-t(1);
for ti = 1:length(t)
    if ti == 1
        du(:,ti) = (A*u0 + I(:,ti) + B*randn(size(B,2),1))*dt;
        u(:,ti) = u0;
    else
        du(:,ti) = (A*u(:,ti-1) + I(:,ti) + B*randn(size(B,2),1))*dt;
        u(:,ti) = u(:,ti-1) + du(:,ti);
    end
end
        