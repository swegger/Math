% Two node simulation package - Brief Tutorial
%
% To be able to use the simulation package, you need to
% first compile the C-programs. In the Matlab command window,
% enter
%
%     mex ccspikesim.cc;
%     mex cciafsim.cc;
%     mex ccsynsim.cc;
%
% ===== Animations similar to Fig. 3
%
% To run a simulation of the two-interval discrimination
% task in the mean-field model with f1=4 and f2=2, type
%
%     mfsim( 4, 2 );
%
% By default, the allowed values of f1 and f2
% range from 1..7. The stimulus frequencies in the main
% text therefore map as follows:
% [10 14 18 22 26 30 34] Hz -> [1 2 3 4 5 6 7]
% Type 'help mfsim' to see more options, e.g., how
% to control the speed of the animation.
%
% ===== Loop through (f1,f2) combinations and create Fig. S4
%
% To run the mean-field simulation for a set of (f1,f2)
% values, just as in Fig. S4, type
%
%     mfmaster;
%
% ===== Loop through (f1,f2) combinations and create Fig. S4
%
% To run a simulation of the spiking network for a set of (f1,f2)
% values, just as in Fig. S4, type
%
%     spikemaster;
%
% Note: by default, 500 Neurons are used per node
% (parameter par.Nneurons in initpar.m) which may take a bit,
% depending on the computational power of your computer.
% To speed things up, just reduce this number---this will
% make the simulations more noisy, though!
%
% ===== Effects of noise in the simulations
%
% To include noise in the mean-field dynamics, change the var
% 'par.mfnoise' in initpar.m (e.g., set par.mfnoise=0.3)
% and run mfsim or mfmaster as before. In particular, with
% noise turned on, check the randomness of decisions if f1==f2,
% e.g.
%
%     mfsim( 4, 4 );
%
% To change the noise level in the spiking dynamics, change
% the number of neurons per node ('par.Nneurons' in initpar.m)
% and rerun spikemaster.
%
% === stability of memory
%
% change the length of the delay period by setting
%
%     par.TSO2 = 7000;
%     par.TSF2 = 7500;
%     par.T    = 8000;
%
% and rerun mfsim/mfmaster/spikemaster. The stability of
% the synaptic tuning can be investigated by changing the
% inhibitory synaptic weight ('par.wI' in initpar.m).
% For zero noise, the memories are stable within +- 0.25%
% To obtain the increased GABA effect (Fig. 4c), increase
% 'par.wI' by 7% and set the noise level to 'par.mfnoise=0.4'
%
% ===== Tuning the i/o function
%
% The i/o function can be recomputed using 'fcurve.m'
% To change the parameters of the neurons and synapses,
% just change the appropriate values in 'initpar.m'.
% Type 'help fcurve' for more details.
%
% The function 'autotuner.m' allows to compute
% i/o functions for a whole set of noise levels and
% synaptic saturation value. For each combination of
% these parameters, 'autotuner' finds the inhibitory
% synaptic weight that yields the best overlap of
% the i/o functions in maintenance mode.
% Type 'help autotuner' for more details.
%
% ==== Contents
%
% see the file 'Contents.m' for a complete list of
% all files. All matlab-files in this package
% are documented, use the help function for more
% information on any one of them.
%
% ====
% This code is part of the supplementary material
% to the report 'Flexible control of mutual inhibition:
% a neural model of two-interval discrimination' by
% Christian K. Machens, Ranulfo Romo, and Carlos D. Brody
