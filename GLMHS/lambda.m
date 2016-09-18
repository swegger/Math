function l = lambda(z,mu_c,C)
%% lambda
%
%   l = lambda(z,B)
%   
%   Computes the activation function of a set of Poisson spiking neurons
%   with a soft threshold as a function of its input, z, and the linear
%   mixing terms, B.
%
%%

l = exp(mu_c + C*z);