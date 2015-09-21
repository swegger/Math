function z = myZscore(x,vec)
%% myZscore
%
%
%%

% for i = 1:size(x,1)
%     X = cat(1,x{i,:});
%     for j = 1:size(x,2)
%         m = mean(X,1);
%         stdev = std(X,[],1);
%         z{i,j} = (x{i,j}-repmat(m,size(x{i,j},1),1))./repmat(stdev,size(x{i,j},1),1);       % Default to normalize across time
%         z{i,j}(isnan(z{i,j})) = 0;      % TODO: is this right?
%     end
% end
% 
% for i = 1:size(x,1)
%     for j = 1:size(x,2)
%         m = mean(x{i,j}(vec{i,j},:),1);
%         stdev = std(x{i,j}(vec{i,j},:),[],1);
%         z{i,j} = (x{i,j}-repmat(m,size(x{i,j},1),1))./repmat(stdev,size(x{i,j},1),1);       % Default to normalize across time
%         z{i,j}(isnan(z{i,j})) = 0;      % TODO: is this right?
%     end
% end


for i = 1:size(x,1)
    for j = 1:size(x,2)
        z{i,j} = sqrt(x{i,j})-repmat(mean(sqrt(x{i,j}),1),size(x{i,j},1),1);
    end
end