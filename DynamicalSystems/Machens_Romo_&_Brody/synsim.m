%synsim - simulates a synapse
%
%   [s,t] = synsim( par, T, spikes ) simulates a
%   synapse that receives an input spike train
%   'spikes' (an array of spike times in msec)
%   for 'T' msec. The synaptic parameters are
%   passed in the structure 'par' with entries
%
%      par.Isatmax = synaptic saturation
%      par.Itau    = synaptic time constant (ms)
%
%   synsim returns the synaptic variable 's' at
%   the time points 't'.
%
%   Note: This "wrapper" function is not autonomous but
%   needs the compiled version of the C-program ccsynsim

%   (c) 2004 CK Machens & CD Brody

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s,t] = synsim( par, T, spikes )

G.dt =  0.1;          %in msec
G.t   = 0:G.dt:T;
G.Nt  = length(G.t);

if (isempty(spikes))
  s = zeros(1,G.Nt);
  return;
end

%integrate-and-fire neuron parameters
G.satmax      =    par.Isatmax;
G.tau         =    par.Itau;

%input and output variables
G.nspikes  = length(spikes);
G.spikes  = spikes;
G.s       = zeros(1,G.Nt);

ccsynsim(G);
s = G.s;
t = G.t;