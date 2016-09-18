function tp = tp_expectation(ts,wm,wp,N)
%% ta_expectation
%
%   ta = ta_expectation(ts,wm)
%
%   Determines the expected value of the production time across all
%   possible measurements and production noise based on the weber fraction
%   of the measurement, wm, and the weber fraction of the production, wp,
%   for N measurements.
%
%   swe 20140321
%       Notes -
%           TODO: Only works for N = 1:2...
%
%%

if nargin < 3
    error('Not enough input arguments!')
elseif N > 2
    error('Values of N > 2 not yet supported!')
end

if N == 1
    % BLS solution for a given tm
    f = @(m,xmin,xmax,wm)(integral(@(x)(x.*exp(-1./2.*(x-m).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)(exp(-1./2.*(x-m).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true));
    
    % Find the average ta for each ts by integrating f(tm)*p(tm|ts) over all tm
    integrand = @(tp,ts,tm,wm,wp,tsmin,tsmax)( tp .* (1./sqrt(2*pi*wp.^2.*f(tm,tsmin,tsmax,wm).^2)) .* exp(-1./2.*(f(tm,tsmin,tsmax,wm)-tp).^2./(wp.*f(tm,tsmin,tsmax,wm)).^2) .* (1./sqrt(2*pi*wm.^2*ts.^2)) .* exp(-1./2.*(tm-ts).^2./(wm.*ts).^2) );
    for i = 1:length(ts)
        tsi = ts(i);
        tp(i) = integral2(@(tp,tm)(integrand(tp,tsi,tm,wm,wp,ts(1),ts(end))),ts(1)-15*wm*ts(1),ts(end)+15*wm*ts(end),ts(1)-15*wm*ts(1),ts(end)+15*wm*ts(end));
    end
    
elseif N == 2
    % BLS solution for a given combination of tm1 and tm2
    f = @(m1,m2,xmin,xmax,wm)(integral(@(x)(x.*(1./sqrt(2.*pi.*wm.^2.*x.^2)).^2.*exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true)./integral(@(x)((1./sqrt(2.*pi.*wm.^2.*x.^2)).^2.*exp(-1./2.*(x-m1).^2./wm.^2./x.^2).*exp(-1./2.*(x-m2).^2./wm.^2./x.^2)),xmin,xmax,'ArrayValued',true));
    
    % Find the average ta for each ts by integrating f(tm1,tm2)*p(tm1,tm2|ts) over all tm1 and tm2
    integrand = @(tp,ts,tm1,tm2,wm,wp,tsmin,tsmax)( tp .* (1./sqrt(2*pi*wp.^2.*f(tm1,tm2,tsmin,tsmax,wm).^2)) .* exp(-1./2.*(f(tm1,tm2,tsmin,tsmax,wm)-tp).^2./(wp.*f(tm1,tm2,tsmin,tsmax,wm)).^2) .* (1./sqrt(2.*pi.*wm.^2.*ts.^2)).^2.*exp(-1./2.*(ts-tm1).^2./wm.^2./ts.^2).*exp(-1./2.*(ts-tm2).^2./wm.^2./ts.^2) );
    for i = 1:length(ts)
        tsi = ts(i);
        tp(i) = integral3(@(tp,tm1,tm2)(integrand(tp,tsi,tm1,tm2,wm,wp,ts(1),ts(end))),ts(1)-15*wm*ts(1),ts(end)+15*wm*ts(end),ts(1)-15*wm*ts(1),ts(end)+15*wm*ts(end),ts(1)-15*wm*ts(1),ts(end)+15*wm*ts(end));
    end
end