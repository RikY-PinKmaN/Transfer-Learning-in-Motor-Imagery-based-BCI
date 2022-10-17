% Feature selection from W, CSP matrix
% ChRanking shows the channels rank
% numchannels is number of the selected channels

function [producedfeatur,selectedw,ind]=ftrain(finaldataset,ChRanking,numchannel,nselectedw)
for i=1:numchannel
    a.x(:,i,:)=finaldataset.x(:,ChRanking(1,i),:);
end

a.y=finaldataset.y;
W=csp(a);
if numchannel>6
    selectedw1=[W(1,:);W(2,:);W(3,:);W(end-2,:);W(end-1,:);W(end,:)];
else
    selectedw1=[W(1,:);W(end,:)];
end
ntrial=size(finaldataset.x,3);
selectedw=[selectedw1;nselectedw];
for trial=1:ntrial
    selectedZ=selectedw*(a.x(:,:,trial))';
    % producedfeature is features*number of channel*number of trials
    producedfeatur.x(:,trial)=log(var( selectedZ(:,:)')./sum(var( selectedZ(:,:)')));
end
producedfeatur.y=finaldataset.y;
% Correlation calculation
for i=1:height(producedfeatur.x)
    coef(i,1)=corr2(producedfeatur.x(i,:)',producedfeatur.y);
end
% Choosing top features
[~,ind]=sort(abs(coef),'descend');
Filter=sort(ind(1:6,1));
producedfeatur.x=producedfeatur.x(Filter(:,1),:);
% selectedw=selectedw(ind(1:6),:);
end