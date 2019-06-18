%% TestFit_Lisberger1999
p = [10,2,100,5,15];

f = @(t,a,b,c,d)( d*(t-a).^b ./ ((t-a).^b + c^b));

wt = 0.01;
t = 0:250;
runNum = 100;

%% Simulate Lisberger1999
for runi = 1:runNum
    v(:,1) = f(t,p(1),p(2),p(3),p(4)) + wt*(t-p(1)).*randn(size(t));
    v(:,2) = f(t,p(1),p(2),p(3),p(5)) + wt*(t-p(1)).*randn(size(t));
    
    v(t<p(1),:) = 0;
    
    % Fit
    [P, feval, exitflg, output] = Fit_Lisberger1999(t,v);
    
    Ps(runi,:) = P;
end


%% Plot
figure('Name','Fit performance')
subplot(1,2,1)
plot(t,v(:,1),'ko')
hold on
plot(t,v(:,2),'bo')

fnc(:,1) = f(t,P(1),P(2),P(3),P(4));
fnc(:,2) = f(t,P(1),P(2),P(3),P(5));
fnc(t<P(1),:) = 0;
plot(t,fnc(:,1),'k')
plot(t,fnc(:,2),'b')

subplot(1,2,2)
relDiff = (repmat(p,[size(Ps,1),1])-Ps)./repmat(p,[size(Ps,1),1]);
errorbar(mean(relDiff,1),std(relDiff,[],1),'o')
hold on
plotHorizontal(0);


