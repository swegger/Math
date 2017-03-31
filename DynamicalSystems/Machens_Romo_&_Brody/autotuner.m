%autotuner - estimate optimal noise and saturation
%
%   res = autotuner( T )
%   computes i/o function for a whole set of values
%   for the noise and synaptic saturation parameters,
%   and finds the optimal synaptic weight for every
%   combination. The optimal, fine-tuned synaptic
%   weight is defined as the one leading to the
%   closest overlap of the tuning curves in the
%   maintenance mode.
%
%   The parameter 'T' determines the precision with
%   which i/o functions are computed---it's the
%   time in seconds used to average over the synaptic
%   output of a neuron. Reasonable values include
%   T = 50...1000 seconds. Smaller values of T can
%   lead to noisy, non-monotonic tuning curves in
%   which case the autotuning breaks down.
%
%   autotuner therefore eats lots of computer time!
%
%   autotuner returns the structure 'res' with the
%   following entries:
%
%      res.gaussnoise = optimal neural noise
%      res.satmax     = optimal synaptic saturation
%      res.wI         = optimal synaptic weight
%      res.satmax_grd = vector of tested saturations;
%      res.noise_grd  = vector of tested noises;
%      res.dist_mat   = matrix of distances
%                       with rows=noises,
%                       and cols=saturations.
%      res.wI_mat     = matrix of optimal weights
%                       with rows=noises,
%                       and cols=saturations.
%      res.gI         = input conductance vector
%      res.sout       = cell matrix of synaptic outputs
%                       with rows=noises,
%                       and cols=saturations.

%   (c) 2004 by CK Machens & CD Brody

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = autotuner( T )

initpar;
figure(1);clf; set( gcf, 'Menubar', 'none' );

%--------------parameters
T = T*1000;                %simulation time in msec
%noise_grd = 0.4:0.02:0.8; %set of noise values
%satmax_grd = 1:0.5:16;    %set of saturation values
noise_grd = 0.3:0.1:0.7;   %set of noise values
satmax_grd = 2:2:12;       %set of saturation values
gI = 0:0.00008:0.006;      %set of inhibit. inputs
gE = par.gE;               %excitatory input

Nn = length(noise_grd);
Ns = length(satmax_grd);
opt_wI = zeros(Nn, Ns);
dist   = zeros(Nn, Ns);
sout   = cell(Nn,Ns);

%---------------loop through noise/satmax values
tic;
for i=1:Nn

  %simulate iaf for a given noise
  par.gaussnoise = noise_grd(i);
  for k=1:length(gI)
    [V,spikes{k}] = iafsim( par, T, gE, gI(k) );
  end
  
  %loop over all synaptic saturation values

  for j=1:Ns

    %simulate synapses for a given saturation
    par.Isatmax = satmax_grd(j);
    for k=1:length(gI)
      stmp = synsim( par, T, spikes{k} );
      s(k) = mean(stmp);
    end
    
    %guess a reasonable weight for nullcline overlap
    smax = max(s);
    auto_wI = gI(min(find( s <= 0.01*smax ))) / smax;
    if isempty(auto_wI) auto_wI=1; end;
    
    %find the optimal weight for nullcline overlap
    tdist = zeros(1,1001);
    for k=5000:15000

      %try this weight
      wI = auto_wI * k/10000;

      %limit curve to invertible, monotonically decreasing part
      u = min(find(diff(s)>=0)); 
      if isempty(u), u=length(s); end;
      if (u>=10)
	is = interp1(wI * s(1:u), gI(1:u), gI);
	
	%compute differences
	uu = s(find(~isnan(is)))*wI - is(find(~isnan(is)));
	uu = sort( abs(uu) );
	pt = uu(1:round(0.65*length(uu)));
	tdist(k-4999) = std(pt);
      else
	tdist(k-4999) = Inf;
      end
    end
    [mindist, ind] = min(tdist);
    dist(i,j) = mindist;
    opt_wI(i,j) = auto_wI*(ind+4999)/10000;
    sout{i}{j} = s;
    
    fprintf(1, 'Did i=%d/%d, j=%d/%d, elapsed time=%.2f\n', ...
	    i, Nn, j, Ns, toc);
    
    %print optimal curve for present values of noise/satmax
    subplot(1,2,1); cla;
    plot( 1000*gI, 1000*opt_wI(i,j)*s, ...
	  1000*opt_wI(i,j)*s, 1000*gI );
    xlabel('nS (plus neuron)');
    ylabel('nS (minus neuron)');
    title(sprintf('optimal overlap for noise=%.2f, satmax=%.1f',...
		  par.gaussnoise, par.Isatmax ));
    drawnow;
  end
end

%--------------plot nullcline distances as function
%              of noises/satmaxes
subplot(1,2,2);
pcolor( satmax_grd, noise_grd, dist ); shading flat;
xlabel( 'satmax' );
ylabel( 'noise' );
colorbar;

%--------------find optimal noise/satmax and
%              create output structure
[mincoldist, ind1] = min( dist );
[mindist, ind2] = min( mincoldist );
i = ind1(ind2);
j = ind2;

res.wI         = opt_wI( i,j );
res.gaussnoise = noise_grd( i );
res.Isatmax    = satmax_grd( j );
res.satmax_grd = satmax_grd;
res.noise_grd  = noise_grd;
res.wI_mat     = opt_wI;
res.dist_mat   = dist;
res.sout       = cell(Nn,Ns);

%--------------plot optimal overlap
subplot(1,2,1); cla;
plot( 1000*gI, 1000*res.wI*sout{i}{j}, ...
      1000*res.wI*sout{i}{j}, 1000*gI );
xlabel('nS (plus neuron)');
ylabel('nS (minus neuron)');
title(sprintf('optimal overlap for noise=%.2f, satmax=%.1f',...
	      res.gaussnoise, res.Isatmax ));

