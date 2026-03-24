function acc = WS_bottlneck_cluster(G)
%acc  = WS_bottleneck_cluster(G)
%    
% Performs the topological clustering using the bottleneck distance. The
% method uses the built-in Ward clustering.
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

% Compute matrix whose entries are pairwise Bottleneck distance. 
% The method is explained in Simulation study 2 in 
%
% [1] Generate random modular network introduced in Songdechakraiwut, T. Chung, 
% M.K. 2022 Topological learning for brain networks, Annals of Applied Statistics arXiv: 2012.00675.
% 
% [2] Songdechakraiwut, T., Shen, L., Chung, M.K. 2021 Topological learning and 
% its application to multimodal brain network integration, Medical Image 
% Computing and Computer Assisted Intervention (MICCAI), LNCS 12902:166-176 
%
% INPUT
% con_i: collection of networks. The size of X matrix should be p x p x n_Group
%                            
% OUTPUT
% lossMtx : pairwise distance matrices
%
%
% (C) 2023 Moo K. Chung
%     University of Wisconsin-Madison
% mkchung@wisc.edu 
%
%  Update history
%     2022 April 3

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
%nBirths x nNetworks and nDeaths x nNetworks in each cluster.

Mb=[];
for j=1:nC
    for i=1:nG
        Mb = [Mb b{i,j}(:,3)];
    end
end

Md=[];
for j=1:nC
    for i=1:nG
        Md = [Md d{i,j}(:,3)];
    end
end

bottle0 = pdist2(Mb', Mb', 'chebychev');
bottle1 = pdist2(Md', Md', 'chebychev');


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


ZWS = linkage(bottle0','ward');
cWS = cluster(ZWS, 'Maxclust',nC);
acc.zero = clustering_accuracy(ytrue, cWS);

ZWS = linkage(bottle1','ward');
cWS = cluster(ZWS, 'Maxclust',nC);
acc.one = clustering_accuracy(ytrue, cWS);





