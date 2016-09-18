function [wmFit, wpFit] = FitModel( ...
    responseT, ...
    invSs, ...
    invSpSupport, ...
    targetDist)
%This function fits parameter values for an ideal observer tracking of a
%ball through a constant trajectory path, leading to an interception when
%the ball reaches a designated point. To simplify computations, and to keep
%the relationship between speed and timing distributions linear, I'm using
%inverse speed rather than speed. I'm also assuming that there is no
%additional noise between the time interval measurement and the calculation
%of inverse speed. In this case, we can model the observer as directly
%taking measurements of inverse speed which are subject to the same scalar
%variance parameter (wm) as the time interval measurements. So it
%doesn't matter what the intercue distances are. As of now, I'm leaving
%them in the code though, in case I add an observer for which it matters.

%This fitting program assumes a uniform prior on inverse speed.

%If you want to use this model to fit an RSG task, just pass in 1 for
%'targetDist'.

%responsT (response times) and invSs (inverse velocity values) must be the
%same size.

initGuess = rand(1,2);

params = ...
    fminsearch(@(x) -1*LogPtp_InvSs(x,targetDist,invSs(:),responseT(:),invSpSupport), ...
    initGuess, ...
    optimset('Display','Iter'));
wmFit = params(1);
wpFit = params(2);
end

%% Functions

function logP = LogPtp_InvSs(x,targetDist,invSs,timeProduced,invSpSupport)
%% The probability of a produced interval given a particular sample
wm = x(1);
wp = x(2);
%%
intBounds = [invSpSupport(1) - 15*wm*invSpSupport(1) ...
    invSpSupport(2) + 15*wm*invSpSupport(2)];
tpIntegrand = ...
    @(invSm,tp) Ptp_vInvEstWp(tp,invSm,wp,wm,targetDist,invSpSupport).* ...
    PinvSm_invSsWm(invSm,invSs,wm);
p = integral(@(invSm) tpIntegrand(invSm,timeProduced), ...
    intBounds(1), ...
    intBounds(2), ...
    'ArrayValued', true);
%%
logP = sum(log(p));
end

function invSpEst = InvSpEst(invSmeasure,wm,invSpSupport)
%% Estimate given a particular measurement
posteriorNumerator = ...
    @(invSs,invSm) 1./sqrt(2*pi*(wm*invSs).^2).* ...
    exp(-(invSs - invSm).^2./(2*(wm*invSs).^2));
alpha = 1/integral(@(s) posteriorNumerator(s,invSmeasure), ...
    invSpSupport(1), ...
    invSpSupport(2), ...
    'ArrayValued', true);
teIntegrand = @(invSs,invSm) alpha*invSs./sqrt(2*pi*(wm*invSs).^2).* ...
    exp(-(invSs - invSm).^2./(2*(wm*invSs).^2));
invSpEst = integral(@(s) teIntegrand(s,invSmeasure), ...
    invSpSupport(1), ...
    invSpSupport(2));

end

function p = Ptp_vInvEstWp(tp,invSm,wp,wm,targetDist,invSpSupport)
%% Probability of producing interval given estimate
invSpEst = InvSpEst(invSm,wm,invSpSupport);
p = 1./sqrt(2*pi*(wp*targetDist*invSpEst).^2).* ...
    exp(-(targetDist*invSpEst - tp).^2./(2*(wp*targetDist*invSpEst).^2));

end

function p = PinvSm_invSsWm(invSm,invSs,wm)
%% Probability of a given measurement given the actual sample
p = 1./sqrt(2*pi*(wm*invSs).^2).*exp(-(invSs - invSm).^2./(2*(wm*invSs).^2));

end
