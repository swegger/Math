%% TestIsomap
%
%
%
%%

% Variables
N = 1000;
neighborRule.type = 'K-closest';
neighborRule.K = 7;

% Generate Data
y = 2*pi*rand(N,2);      % Randomly generates points on a plane spanning 0 to 2*pi
X = [y(:,1).*cos(y(:,1)) y(:,2) y(:,1).*sin(y(:,1))];           % Converts data into a "Swiss roll" in 3-D space

% Compute isomap
Y = Isomap(X,2,'neighborRule',neighborRule,'TestY',y);        % Computes the lower dimensional coordinates

% Plot the results
figure
plot(y(:,1),y(:,2),'.');
hold on
plot(Y(:,1),Y(:,2),'o');
