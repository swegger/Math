function pieceWiseExpFallFitterTest(varargin)
%% pieceWiseExpRiseFitterTest
%
%   pieceWiseExpRiseFitterTest()
%   
%       Tests identifiability of peice-wise exponential model to data w/
%       Gaussian noise.
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'t',0:500)
addParameter(Parser,'p',[20,0,100,200])
addParameter(Parser,'sigma',1)
addParameter(Parser,'fitN',100)
addParameter(Parser,'p0',[20*randn,10*randn,200*rand,300*rand])
addParameter(Parser,'pNames',{'R_0','R_f','\tau','t_0'})

parse(Parser,varargin{:})

t = Parser.Results.t;
p = Parser.Results.p;
sigma = Parser.Results.sigma;
fitN = Parser.Results.fitN;
p0 = Parser.Results.p0;
pNames = Parser.Results.pNames;

%% Generate data set
T = repmat(t(:),[1,fitN]);
R = repmat(pieceWiseExpFall(t(:),p),[1,fitN]) + sigma*randn(size(T));

%% Fit
for fiti = 1:fitN
    phat(fiti,:) = pieceWiseExpFallFitter(t(:),R(:,fiti),p0,[Inf,Inf,Inf,Inf],[-Inf,-Inf,0,0]);
    Rhat(:,fiti) = pieceWiseExpFall(t(:),phat(fiti,:));
end

%% Plot results
figure('Name','Fit parameters verses actual data')
for pi = 1:length(p)
    subplot(2,length(p),pi)
    if p(pi) == 0
        histogram((phat(:,pi)-p(pi)),ceil(0.2*fitN));
        xlabel('Absolute error')
    else
        histogram((phat(:,pi)-p(pi))/p(pi),ceil(0.2*fitN));
        xlabel('Relative error')
    end
    hold on
    plotVertical(0);
    title(pNames{pi})
end

subplot(2,length(p),[length(p)+1 2*length(p)])
plot(t,R,'o','Color',[0.6 0.6 0.6])
hold on
plot(t,Rhat,'k')