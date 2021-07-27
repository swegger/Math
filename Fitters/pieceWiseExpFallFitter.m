function p = pieceWiseExpFallFitter(tin,R,p0,ub,lb)

% p = fminsearch(@(p)minimizer(tstar(:),tref(:),p),p0);

p = fmincon(@(p)minimizant(tin(:),R(:),p),p0,[],[],[],[],...
    lb,ub);

%% Functions
function SumSqErr = minimizant(tin,R,p)
%%
SumSqErr = sum( (pieceWiseExpRise(tin,p) - R).^2 );
    
    
        