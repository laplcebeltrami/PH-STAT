function acc = GH_cluster(G)
%acc =GH_cluster(G)
%    
% Performs hierarchical_cluster(G) on graphs G
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
%
% (C) 2020 Moo K. Chung
%     University of Wisconsin-Madison
%  Contact mkchung@wisc.edu for support 
%
% Update history
%   December 22, 2021 create
%   Otober 24, 2022 updated
%   March 25, 2023 made into a function


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

ZWS = linkage(C','ward');
cWS = cluster(ZWS, 'Maxclust',nC);
acc = clustering_accuracy(ytrue, cWS);

%clustergram(ZWS) %, 'RowLabels', [], 'ColumnLabels', [], 'Linkage', 'complete', 'Dendrogram', 'row')
