function [du, u_null] = NL_linear2DVectorField(u,A,F,I)
%% FitzHughNagumoModelVectorField
%
%
%%

% Vector field
for i = 1:2
    fU(:,i) = F{i}(u);
end
du = A*fU' + I;


% Nullclines
u_null = NaN;       % TODO: does matlab have a solver?
