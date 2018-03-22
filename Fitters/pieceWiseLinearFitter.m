function p = pieceWiseLinearFitter(tin,tout,p0,ub,lb)

% p = fminsearch(@(p)minimizer(tstar(:),tref(:),p),p0);

p = fmincon(@(p)minimizant(tin(:),tout(:),p),p0,[],[],[],[],...
    lb,ub);

%% Functions
function SumSqErr = minimizant(tin,tout,p)
%%
SumSqErr = sum( (pieceWiseLinear(tin,p) - tout).^2 );
    
    
        