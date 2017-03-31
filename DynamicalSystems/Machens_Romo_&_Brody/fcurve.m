%fcurve - computes a neuron's i/o function
%
%   [x, fx, lambda] = fcurve( T ) computes the
%   integrate-and-fire neuron's i/o function,
%   i.e., its mean synaptic output as a function
%   of constant conductance inputs. The parameter
%   'T' specifies the simulation length (in sec)
%   used to estimate the synaptic output. The
%   i/o function is returned in the vectors 'x'
%   and 'fx'; look at it with 'plot( x,fx )'.
%   The returned variable 'lambda' is the
%   excitation/inhibition factor as used in Eq. ???
%   in the supplementary materials.
%
%   fcurve uses the parameters specified in
%   initpar.m. To recompute the i/o function
%   with the same precision as used in the paper,
%   set T=1000 sec.

%   (c) 2004 CK Machens & CD Brody

function [xgrd, fx, lambda] = fcurve( T )

initpar;
T = T*1000;  %convert to milliseconds

%--------------inhibitory conductance inputs
gI  = 0:0.0001:0.006;
figure(1);clf;
set( gcf, 'Menubar', 'none' );
subplot(2,2,1); hold on;
axis( [ gI(1) gI(end) 0 0.75*par.Isatmax ] );

%--------------compute i/o-function for maintenance mode
for i=1:length(gI)
  [V,spikes] = iafsim( par, T, par.gE, gI(i) );
  stmp = synsim( par, T, spikes );
  s(i) = mean( stmp );
  pl1 = plot( gI(1:i), s(1:i),  'k.-' );
  drawnow;
end

%--------------compute i/o function for loading mode
for i=1:length(gI)
  [V,spikes] = iafsim( par, T, par.gE_load, gI(i) );
  stmp = synsim( par, T, spikes );
  s2(i) = mean( stmp );
  pl2 = plot( gI(1:i), s2(1:i),  'r.-' );
  drawnow;
end

%--------------compute i/o function for decision mode
for i=1:length(gI)
  [V,spikes] = iafsim( par, T, par.gE_deci, gI(i) );
  stmp = synsim( par, T, spikes );
  s3(i) = mean( stmp );
  pl3 = plot( gI(1:i), s3(1:i),  'b.-' );
  drawnow;
end

%--------------finish plot
legend( [pl2, pl1, pl3], 'Loading', 'Maintenance', 'Decision' );
xlabel('inhibitory conductance input (uS)');
ylabel('synaptic output');
title('i/o functions for three excitatory inputs')

%--------------show maintenance mode overlap
subplot(2,2,2);
plot( gI, par.wI*s, par.wI*s, gI );
xlabel('uS (plus neuron)');
ylabel('uS (minus neuron)');
title('memory mode for present synaptic weights')
ax = axis;
axis( [-0.0005 ax(2) -0.0005 ax(4)] );

%--------------estimate lambda, loading->maintenace mode
gmove = 0:0.00005:0.004;
for k=1:length(gmove)
  ind = min(find( s <= 0.01*max(s) ));
  ss2 = interp1( gI-gmove(k), s2, gI(1:ind));
  ds(k) = std( s(1:ind)-ss2 );
end
[dsmin, dsind] = min( ds );
dgE = par.gE_load-par.gE;
lambda = gmove(dsind) / dgE;

subplot(2,2,3);
pl1 = plot( gI, s, 'k' ); hold on;
pl2 = plot( gI-dgE*lambda, s2, 'r' );
dgE = par.gE_deci-par.gE;
pl3 = plot( gI-dgE*lambda, s3, 'b' );
axis( [gI(1)-gmove(end) gI(end) 0 max(s2)] );
legend( [pl2, pl1, pl3], 'Loading', 'Maintenance', 'Decision', 3 );
xlabel('inhibitory conductance input (uS)');
ylabel('synaptic output');
title('shifted i/o functions (lambda-approximation)')

%--------------flip the curve ...
ggrd = 0:0.00001:0.01;
fx = interp1( gI, s, ggrd );
xgrd = lambda*par.gE - ggrd;
xgrd = xgrd(end:-1:1);
fx   = fx(end:-1:1);

%--------------...and extrapolate the right hand side
fx2 = interp1( gI, s2, ggrd );
ind = find( fx2 > max(fx) );
fx_ext = fx2( ind );
dxgrd = xgrd(2)-xgrd(1);
xgrd = [xgrd, xgrd(end) + dxgrd + xgrd(1:length(ind)) - xgrd(1)];
fx = [fx, fx_ext(end:-1:1)];

%--------------plot the finished i/o function
subplot(2,2,4);
plot( xgrd, fx );
xlabel('conductance input (uS)');
ylabel('synaptic output');
title('generic i/o function')
axis( [min(xgrd) max(xgrd) 0 max(fx)*1.1]);

