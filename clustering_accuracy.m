function [accuracy C]=cluster_accuracy(ytrue,ypred)
% function accuracy=cluster_accuracy(ytrue,ypred)
%
% Find the clustering accuracy of prediction, given the true labels
%
% INPUT
% ytrue = a vector of true labels
% ypred = a vector of the predicted labels

% OUTPUT
% acc = Accuracy of clustering results
%
%
% The function is downloaded from
% https://github.com/laplcebeltrami/PH-STAT
%
% (C) 2021 Moo K. Chung
%     University of Wisconsin-Madison
%  Contact mkchung@wisc.edu for support 
%
% Update history
%   2021 Nov 28 
%   2023 Apri 11 function name changed from clustering_accuracy.m to cluster_accuracy.m


%test data
%ytrue = [ 1 1 1  2 2 3 3]
%ypred = [ 1 1 2  1 1 3 3]

n=length(ytrue); % number of samples

%confusion matrix
C = confusionmat(ytrue,ypred); 

%figure; heatmap(C); colorbar

%The iteration below also computes the confusion matrix
%k=max([ytrue(:); ypred(:)]) %find # of cluster
%C=zeros(k);
%for i=1:n  
%    C(ypred(i),ytrue(i))=C(ypred(i),ytrue(i))+1;
%end

%test example
%2     2     0
%1     0     0
%0     0     2
     

%solve the linear assignment problem of maximizing the clusting result
% 0 is the cost of mismatch. We ony care matched pairs
M=matchpairs(C, 0, 'max'); 

% test example
% 2 1 
% 1 2
% 3 3

accuracy=sum(C(sub2ind(size(C), M(:,1), M(:,2))))/n;

%permuated confusion matrix that maximizes the diagonal elements

% P=[]; %permutation matrix
% for i=1:size(M,1)
%   P(M(i,1),M(i,2)) = 1;
% end
% 
% C=P*C*P;
% figure;  heatmap(C); colorbar

