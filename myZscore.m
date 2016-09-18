function z = myZscore(x,vec)
%% myZscore
%
%
%%

%% Transform according to mean and std across all the data
% for i = 1:size(x,1)
%     X = cat(1,x{i,:});
%     m = mean(X,1);
%     stdev = std(X,[],1);
%     for j = 1:size(x,2)
%         z{i,j} = (x{i,j}-repmat(m,size(x{i,j},1),1))./repmat(stdev,size(x{i,j},1),1);       % Default to normalize across time
%         z{i,j}(isnan(z{i,j})) = 0;      % TODO: is this right?
%     end
% end


%% Transform according to values of x specified by vector input
% for i = 1:size(x,1)
%     for j = 1:size(x,2)
%         m = mean(x{i,j}(vec{i,j},:),1);
%         stdev = std(x{i,j}(vec{i,j},:),[],1);
%         z{i,j} = (x{i,j}-repmat(m,size(x{i,j},1),1))./repmat(stdev,size(x{i,j},1),1);       % Default to normalize across time
%         z{i,j}(isnan(z{i,j})) = 0;      % TODO: is this right?
%     end
% end

%% Square-root transform
% for i = 1:size(x,1)
%     for j = 1:size(x,2)
%         z{i,j} = sqrt(x{i,j})-repmat(mean(sqrt(x{i,j}),1),size(x{i,j},1),1);
%     end
% end

%% Normalize to peak firing rate across conditions
X = cat(1,x{:,:});
maxX = max(X,[],1);
maxX(maxX == 0) = 1;
for i = 1:size(x,1)
    for j = 1:size(x,2);
        z{i,j} = x{i,j}./repmat(maxX,size(x{i,j},1),1);
    end
end