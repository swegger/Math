function [du] = linear2DVectorField(u,A,I)
%% FitzHughNagumoModelVectorField
%
%
%%

% Vector field
for i = 1:size(u,1)
    for j = 1:size(u,2)
        du(i,j,:) = A*permute(u(i,j,:),[3,2,1]) + I;
    end
end