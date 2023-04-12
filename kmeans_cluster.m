function acc = kmeans_cluster(G)
%function acc = kmeans_cluster(G)
%    
% Performs kmeans clustering on a collection of graphs G
%
% INPUT
% G : collection of graphs stored as the cell of size # of graphs x # of clusters
%
% OUTPUT
% acc: clustering accuracy
%
% The topologial cluserting method is explained in 
% 
% [1] % Chung, M.K., Huang, S.-G., Carroll, I.C., Calhoun, V.D., Goldsmith, H.H. 
% 2023 Topological  State-Space Estimation of Functional Human Brain Networks. arXiv:2201:00087
%
% [2] Moo K. Chung, Camille Garcia Ramos, Felipe Branco De Paiva, Jedidiah Mathis, Vivek Prabharakaren, 
% Veena A. Nair, Elizabeth Meyerand, Bruce P. Hermann, Jeffery R. Binder, Aaron F. Struck, 2023 
% Unified Topological Inference for Brain Networks in Temporal Lobe Epilepsy Using the Wasserstein Distance, 
% arXiv:2302.06673
%
% 
% The topologial distace is given in 
%
% [3] Songdechakraiwut, T., Shen, L., Chung, M.K. 2021 Topological learning and 
% its application to multimodal brain network integration, Medical Image 
% Computing and Computer Assisted Intervention (MICCAI), LNCS 12902:166-176 
%
% [4] Songdechakraiwut, T. Chung, M.K. 2023 Topological learning 
% for brain networks, Annals of Applied Statistics 17:403-433, arXiv: 2012.00675
%
% If you are using any part of the code, please reference the above paper.
% The function is downloaded from 
% http://pages.stat.wisc.edu/~mchung/publication.html
%
% The function is downloaded from
% https://github.com/laplcebeltrami/PH-STAT
%
%
% (C) 2020 Moo K. Chung
%     University of Wisconsin-Madison
%  Contact mkchung@wisc.edu for support 
%
% Update history
%   December 22, 2021 create
%   Otober 24, 2022 updated
%   March 25, 2023 made into a function


%% Compute set of births and set of deaths

nG = size(G,1); %number of graphs in each cluster
nC = size(G,2); %number of clusters

nNetworks = nG*nC; %total number of networks

nNodes = size(G{1,1},1); % # of nodes
nEdges = nNodes*(nNodes-1)/2; % # of edges

%C=zeros(nEdges,nNetworks); 
C=[];
for j=1:nC
    for i=1:nG
        edgeweights = adj2vec(G{i,j});
        C=[C  edgeweights'];
    end
end

ytrue=[]; %true labels
for i = 1:nC
    % add i to the sequence of clusters
    ytrue = [ytrue; repmat(i, nG, 1)];
end

accuracy = [];
for i=1:100
    ypred = kmeans(C',nC);
    [clustering_acc Conf] = clustering_accuracy(ytrue,ypred);
    accuracy(i)=clustering_acc;
end

acc.mean= mean(accuracy);
acc.std = std(accuracy);


%------------------------
function vec = adj2vec(adj);
%
% function vec = adj2vec(adj);
% vectorize the square matrix and produce the vector of elemenets 
% in upper triangle. For example, from the distance matrix, it produces the pdist.
% 
% INPUT:
% adj has to be bigger than 2 x 2 matrix
%
%
% (C) 2019 Moo K. Chung 
% University of Wisconsin-Madison
% mkchung@wics.edu
%

n=size(adj,1);
vec=adj(1,2:end);
 
if n>=2
    for i=2:n
        vec= [vec adj(i,i+1:end)];
    end;
    
end


%-------------
function [accuracy C]=clustering_accuracy(ytrue,ypred)
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
% (C) 2021 Moo K. Chung
%     University of Wisconsin-Madison
%  Contact mkchung@wisc.edu for support 
%
% Update history
%   2021 Nov 28 
%   2023 Apri 11 function name changed from clustering_accuracy.m to cluster_accuracy.m

n=length(ytrue); % number of samples
C = confusionmat(ytrue,ypred); 
M=matchpairs(C, 0, 'max'); 
accuracy=sum(C(sub2ind(size(C), M(:,1), M(:,2))))/n;
