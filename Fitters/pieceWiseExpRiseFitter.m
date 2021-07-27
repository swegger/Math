function p = pieceWiseExpRiseFitter(tin,R,p0,ub,lb)

% p = fminsearch(@(p)minimizer(tstar(:),tref(:),p),p0);

options = optimset('Display','off');
p = fmincon(@(p)minimizant(tin(:),R(:),p),p0,[],[],[],[],...
    lb,ub,[],options);

%% Functions
function SumSqErr = minimizant(tin,R,p)
%%
SumSqErr = sum( (pieceWiseExpRise(tin,p) - R).^2 );
    
    
        