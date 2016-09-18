function [mu, V] = kalmanPredictor(mu0,V0,A,Gamma,K)
%% kalmanPredictor
%
%   [mu, V] = kalmanPredictor(mu0,V0,A,Gamma,K)
%
%   Finds the mean, mu, and variance, V, of the predictive distribution of
%   latentent variable, z, for K steps from intitial conditions mu0, V0 and
%   using the state equation
%       z(k+1) = Az(k) + eta
%   Where eta is a random varible distributed as a Gaussian with covariance
%   Gamma.
%
%%

for k = 1:K
    if k == 1
        mu(:,k) = A*mu0;
        V(:,:,k) = A*V0*A' + Gamma;
    else
        mu(:,k) = A*mu(:,k-1);
        V(:,:,k) = A*V(:,:,k-1)*A' + Gamma;
    end
end