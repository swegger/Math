function [mu, r, J] = KmeansClustering(x,K,mu_init,tolerance)
%% KmeanClustering
%
%   [mu, r] = KmeansClustering(x,K,mu_init,r_init);
%
%   Assigns data points in x to one of K clusters that minimize the mean
%   squared distance of each point in a cluster from that cluster's mean.
%
%%

% Assign each data point to a cluster
mu = mu_init;
for i = 1:K
    js(:,i) = sum((x - repmat(mu(i,:),size(x,1),1)).^2,2);
end
J = sum(js(:));
Jold = J*2;

disp(J)
while (Jold-J) > tolerance
    Jold = J;
    
    % Assign each data point to a cluster
    for i = 1:K
        js(:,i) = sum((x - repmat(mu(i,:),size(x,1),1)).^2,2);
    end
    r = (js-repmat(min(js,[],2),1,size(js,2)) == 0);
    
    % Find the new mean of each cluster
    for i = 1:K
        mu(i,:) = sum(repmat(r(:,i),1,size(x,2)).*x,1)/sum(r(:,i));
    end
    
    % Determine the new distortion measure
    for i = 1:K
        js(:,i) = sum((x-repmat(mu(i,:),size(x,1),1)).^2,2);
    end
    J = sum(js(:));
    disp(J);
end
        
