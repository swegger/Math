function [u, du] = NL_linear2D(t,u0,A,F,I,B)
%% linear2D
%%

if length(I) == 1
    I = I*ones(size(t));
end

dt = t(2)-t(1);
for ti = 1:length(t)
    if ti == 1
        for i = 1:2
            fU(:,i) = F{i}(u0);
        end
        du(:,ti) = (A*fU+ I(:,ti) + B*randn(size(B,2),1))*dt;
        u(:,ti) = u0;
    else
        for i = 1:2
            fU(:,i) = F{i}(u(:,ti-1));
        end
        du(:,ti) = (A*fU + I(:,ti) + B*randn(size(B,2),1))*dt;
        u(:,ti) = u(:,ti-1) + du(:,ti);
    end
end
        