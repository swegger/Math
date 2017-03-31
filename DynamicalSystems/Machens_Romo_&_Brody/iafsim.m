%iafsim - simulates integrate-and-fire neuron
%
%   [V,spikes,t] = iafsim( par, T, gE, gI ) simulates
%   the integrate-and-fire neuron for 'T' msec and
%   returns the voltage trace 'V' at the time
%   points 't' and the spike times 'spikes'
%   in msec. The parameters of the integrate-and-fire
%   neurons are passed in the structure 'par' with
%   entries
%
%      par.C        = capacitance (nF)
%      par.Vthresh  = treshold (mV)
%      par.Vreset   = reset potential (mV)
%      par.refrac   = refractory period (ms)
%      par.EL       = cell's resting potential (mV)
%      par.EE       = rev. potential of excit. input (mV)
%      par.EI       = rev. potential of inhib. input (mV)
%      par.gL       = leak conductance (nS)
%      par.gaussnoise = std. dev. of additive noise (mV)
%
%   The variables 'gE' and 'gI' denote the constant
%   excitatory and inhibitory conductance input to
%   the neuron.
%
%   Note: This "wrapper" function is not autonomous but
%   calls the compiled version of the C-program cciafsim.cc

%   (c) 2004 by CK Machens & CD Brody

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [V,spikes,t] = iafsim( par, T, gE, gI )

%--------------simulation parameters
G.dt =  0.1;          %in msec
G.t   = 0:G.dt:T;
G.Nt  = length(G.t);

%--------------integrate-and-fire neuron parameters
G.C           =    par.C;
G.Vthresh     =    par.Vthresh;
G.Vreset      =    par.Vreset;
G.refrac      =    par.refrac;
G.EL          =    par.EL;
G.EE          =    par.EE;
G.EI          =    par.EI;
G.gL          =    par.gL;
G.gE          =    gE;
G.gI          =    gI;
taum          =    par.C/par.gL;
var           =    par.gaussnoise.^2 * (1-G.dt/taum)/(G.dt/taum);
G.gaussnoise  =    sqrt(var)*0.01418;

%--------------pregenerated random numbers
G.grand       =    randn(1,G.Nt);

%--------------reporting variables
G.V               = G.EL*ones(1,G.Nt);
G.nspikes         = ceil(T/par.refrac);
G.spikes          = zeros(1,G.nspikes);

cciafsim(G);

if (G.nspikes>=1)
  spikes = G.spikes(:,1:G.nspikes);
else
  spikes = [];
end
V = G.V;
t = G.t;