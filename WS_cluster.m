function acc = WS_cluster(G)
%function acc = WS_cluster(G)
%    
% Performs the topological clustering using the Wassersetin distance. The
% method uses the built-in k-means clustering which depends on intial
% seeds. Thus, need to perform the method multiple times with different
% seeds. 
%
% INPUT
% G : collection of graphs stored as the cell of size # of graphs x # of clusters
%
% OUTPUT
% acc: clustering accuracy based on 100 inpendent clustering reporting mean
% and standard deviation.
%
% The topologial cluserting method is explained in 
% 
% [1] % Chung, M.K., Huang, S.-G., Carroll, I.C., Calhoun, V.D., Goldsmith, H.H. 
% 2023 Topological  State-Space Estimation of Functional Human Brain Networks. arXiv:2201:00087
% 
% The topologial distace is given in 
%
% [2] Songdechakraiwut, T., Shen, L., Chung, M.K. 2021 Topological learning and 
% its application to multimodal brain network integration, Medical Image 
% Computing and Computer Assisted Intervention (MICCAI), LNCS 12902:166-176 
%
% [3] Songdechakraiwut, T. Chung, M.K. 2023 Topological learning 
% for brain networks, Annals of Applied Statistics 17:403-433, arXiv: 2012.00675
%
% If you are using any part of the code, please reference the above paper.
% The function is downloaded from 
% http://pages.stat.wisc.edu/~mchung/publication.html
%
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

b=cell(nG,nC); % # of networks x # of clusters
d=cell(nG,nC);

% Perform birth-death decomposition 
for j=1:nC 
    for i=1:nG
        [b{i,j}, d{i,j}] =WS_decompose(G{i,j});
    end
end

%Tranforms birth and death values into a matrix form of size
%(nBirths+nDeaths) x nNetworks in each cluster.

M=[];
for j=1:nC
    for i=1:nG
        bd=[b{i,j}(:,3); d{i,j}(:,3)];
        M = [M bd];
    end
end

%M should be stored as [cluster1, cluster2, ...]
%M is a matrix of size # of (brith + death values) x # of subjects in each
%cluster

%The clustering_accuracy.m computation is based on the techncial report 
%https://github.com/laplcebeltrami/clustering/blob/main/clustering.accuracy.pdf

ytrue=[]; %true labels
for i = 1:nC
    % add i to the sequence of clusters
    ytrue = [ytrue; repmat(i, nG, 1)];
end

accuracy=[];
for i=1:100
    ypred = kmeans(M',nC);
    [clustering_acc C] = clustering_accuracy(ytrue,ypred);
    accuracy(i)=clustering_acc;
end

acc.mean= mean(accuracy);
acc.std = std(accuracy);


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



%----------------
%OLD CODE

% n1 = length(g1);
% n2 = length(g2);
% n3 = length(g3);
% n4 = length(g4);
% 
% nNetworks = n1+n2+n3+n4; %total number of networks
% 
% % Perform birth-death decomposition 
% for i=1:n1
% [b1{i}, d1{i}] =WS_decompose(g1{i});
% end
% 
% for i=1:n2
% [b2{i}, d2{i}] =WS_decompose(g2{i});
% end
% 
% for i=1:n3
% [b3{i}, d3{i}] =WS_decompose(g3{i});
% end
% 
% for i=1:n4
% [b4{i}, d4{i}] =WS_decompose(g4{i});
% end
% 
% %Tranforms death values into a matrix form of size [ (nBirths+nDeaths) x
% nNetworks, (nBirths+nDeaths) x nNetworks, (nBirths+nDeaths) x
% nNetworks...]
% 
% M1=[];
% for i=1:n1
%     M1(:,i)=[b1{i}(:,3); d1{i}(:,3)];
% end
% 
% M2=[];
% for i=1:n2
%     M2(:,i)=[b2{i}(:,3); d2{i}(:,3)];
% end
% 
% M3=[];
% for i=1:n2
%     M3(:,i)=[b3{i}(:,3); d3{i}(:,3)];
% end
% 
% M4=[]; %cluster 4
% for i=1:n2
%     M4(:,i)=[b4{i}(:,3); d4{i}(:,3)];
% end
%
% M4 is stored as   [subj1 sub2 .....] of size (nBirths+nDeaths) x nNetworks


% M=[M1, M2, M3, M4];
% 
% %M is a matrix of size # of (brith + death values) x # of subjects
% %The clustering_accuracy.m computation is based on the techncial report 
% %https://github.com/laplcebeltrami/clustering/blob/main/clustering.accuracy.pdf
% accuracy=[];
% for i=1:100
%     ypred = kmeans(M',4);
%     ytrue = [ones(5,1); 2*ones(5,1); 3*ones(5,1); 4*ones(5,1)];
%     [acc C] = clustering_accuracy(ytrue,ypred);
%     accuracy(i)=acc;
% end