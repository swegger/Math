% Two node simulation package - Contents
%
% Contents.m    - this file
% Tutorial.m    - brief tutorial. type 'help Tutorial'
%
% -------------Matlab m-files
%
% autotuner.m   - find optimal satmax and noise values
% fcurve.m      - compute neuron's i/o function
% iafsim.m      - simulate integrate-and-fire neuron
% initpar.m     - initialize parameters used in two-node model
% mfmaster.m    - simulate mean-field model for several
%                 (f1,f2)-values and creates Fig. S???
% mfsim.m       - animated simulation of two-node
%                 mean-field network, similar to Fig. 3
% rainbow_colors.m - auxiliary function
% spikemaster.m - simulate spiking model for several
%                 (f1,f2)-values and creates Fig. S???
% spikesim.m    - simulate two-node spiking network
% synsim.m      - simulate synapse
%
% -------------Matlab mat-files
%
% fcurve.mat    - precomputed i/o function 
%
% -------------C-programs
%
% cciafsim.cc   - engine for integrate-and-fire simulation
% ccsynsim.cc   - engine for synapse simulation
% ccspikesim.cc - engine for spiking network simulation
%
% -------------
% This code is part of the supplementary material
% to the report 'Flexible control of mutual inhibition:
% a neural model of two-interval discrimination' by
% Christian K. Machens, Ranulfo Romo, and Carlos D. Brody
