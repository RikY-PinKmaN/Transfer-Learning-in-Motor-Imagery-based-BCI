function [producedfeatur]=ftest(finaldataset,ChRanking,numchannel,selectedw,ind)
for i=1:numchannel
    a.x(:,i,:)=finaldataset.x(:,ChRanking(1,i),:);
end

a.y=finaldataset.y;
ntrial=size(finaldataset.x,3);
for trial=1:ntrial
    selectedZ=selectedw*(a.x(:,:,trial))';
    % producedfeature is features*number of channel*number of trials
    producedfeatur.x(:,trial)=log(var( selectedZ(:,:)')./sum(var( selectedZ(:,:)')));

end
producedfeatur.y=finaldataset.y;
Filter=sort(ind(1:6,1)); % Filtering based on correlation
producedfeatur.x=producedfeatur.x(Filter(:,1),:);