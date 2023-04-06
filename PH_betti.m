function  beta = PH_betti(C, thresholds)
%function beta = PH_betti(C, thresholds)
%
% The function computes the 0-th (the number of connected components) and 1st Betti numbers 
% (cycles) and the size of the largest cycle over the range of filtration 
% values in the 1-skeleton of connectivity matrix C
%
% INPUT
% C           : Weighted connectivity matrix or edge weights
% thresholds  : Range of filtration/thresholds to use in C. 
%
%               Example. [0:0.01:1] will tresholds C between [0,1] at 0.01 increment. 
%
%
% OUTPUT
% beta        : beta.zero is the 0th Betti curve 
%               beta.one is the 1st Betti curve 
%
%          
% The matheamtical details of the methods are published in [1] and [2]. 
% If you are using this code, please reference [1] or [2]. 
%
%% [1] Chung, M.K., Lee, H. Ombao. H., Solo, V. 2019 Exact topological inference 
%%     of the resting-state brain networks in twins, Network Neuroscience 3:674-694
%%     http://www.stat.wisc.edu/~mchung/papers/chung.2019.NN.pdf
%
% [2] Chung, M.K., Huang, S.-G., Gritsenko, A., Shen, L., Lee, H. 2o19
%    Statistical inference on the number of cycles in brain networks. 
%    IEEE International Symposium on Biomedical Imaging (ISBI) 113-116 
%
% Given two Betti curves, we can perform the Exact Topologica Inference (ETI) that provide
% statistical signficance (p-value) of testing the equivalencce of two
% curves. ETI is introduced in
%
% [3] Chung, M.K., Vilalta, V.G., Lee, H., Rathouz, P.J., Lahey, B.B., Zald, D.H. 
%     2017 Exact topological inference for paired brain networks via persistent 
%     homology. Information Processing in Medical Imaging (IPMI) 10265:299-310
%     http://www.stat.wisc.edu/~mchung/papers/chung.2017.IPMI.pdf
%
% [4] Chung, M.K., Luo, Z., Leow, A.D., Alexander, A.L., Richard, D.J., Goldsmith, H.H. 
%     2018 Exact Combinatorial Inference for Brain Images, Medical Image Computing and 
%     Computer Assisted Intervention (MICCAI), 11070:629-637
%     http://www.stat.wisc.edu/~mchung/papers/chung.2018.MICCAI.pdf
%
% [5] Chung, M.K. Lee, H., Gritsenko, A., DiChristofano, A., Pluta, D. 
%     Ombao, H. Solo, V. Topological Brain Network Distances, ArXiv 1809.03878
%     http://arxiv.org/abs/1809.03878
%
%
% (C) 2017- Moo K. Chung, Hyekyung Lee                           
%      University of Wisconsin-Madison
%      Seuoul National University
%      mkchung@wisc.edu
%
% Update history    
% 2017 December 18. errors fixed
% 2018 Jun.  16  nargin added
% 2018 Jul.  28  Betti-1 number computation added. 
% 2018 Aug.  11  con2adj error fixed. code simplified using built-in MATLAB function
% 2023 April 5  outputs simplified. 


if nargin<=1
    %if threshold is not given, thresholds are automatically set between
    %the minimum and maximum of edge weghts at 0.01 increment.
    maxC= max(max(C));
    minC = min(min(C));
    thresholds= minC:0.01:maxC;
end


%-------------------------------
% The algorithm requries computing beta0 first. Then computes beta1 using beta0
% beta1 is a function of beta0. See [1] or [2] for details.
% Note Euler characteristic = beta0 - beta1 = # of nodes - # of edges. 
% Thus, beta1 = beta0 - # of nodes + # of edges


beta0 =[];
biggest0=[];

beta1 =[];
biggest1=[];

n_nodes = size(C,1);

%C=abs(C); %removes possible negative values in the connectivity matrix. 

for rho=thresholds  %this range needs to be changed depending on the value of C
    %computest Beti-0
    %computes adjacency matrix    
    adj = sparse(C>rho); %introduces diagonal entries
    adj = adj - diag(diag(adj)); %removes diagonal entries
    
    [n_components,S] = conncomp(adj); %faster routine
    %[n_components,S] = graphconncomp(adj); built-in MATLAB function is slow
    %n_components is the number of components
    %S is a vector indicating to which component each node belongs
    beta0=[beta0 n_components]; %Betti_0: the number of connected components
    
    nn = hist(S,[1:n_components]); 
    % nn contains the number of nodes in each connected component 
    biggest0 = [biggest0 max(nn)];

    %computes Beti-1
    n_edges=sum(sum(adj))/2;
    n_cycle = n_components - n_nodes + n_edges;
    beta1=[beta1 n_cycle];
    
    %if n_cycle<=0
    %    n_cycle=0; % the number of cycle may go below zero numerically. 
    %end;
    
end

beta.zero =beta0;
beta.one  =beta1;



function [S,C] = conncomp(G)
  % CONNCOMP Drop in replacement for graphconncomp.m from the bioinformatics
  % toobox. G is an n by n adjacency matrix, then this identifies the S
  % connected components C. This is also an order of magnitude faster.
  %
  %  http://www.alecjacobson.com/weblog/?p=4203
  %
  %
  % [S,C] = conncomp(G)
  %
  % Inputs:
  %   G  n by n adjacency matrix
  % Outputs:
  %   S  scalar number of connected components
  %   C  
  
  [p,q,r] = dmperm(G+speye(size(G)));  %Dulmage-Mendelsohn permutation.
  S = numel(r)-1;
  C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
  C(p) = C;


% OLD CODE
% for rho=thresholds  %this range needs to be changed depending on the value of C
%     %computest Beti-0
%     adj = con2adj(C,rho);
%     [n_components,sizes,members] = networkComponents(adj);
%     
%     %computes Beti-1
%     n_edges=sum(sum(adj))/2;
%     n_cycle = n_components - n_nodes + n_edges;
%     if n_cycle<=0
%         n_cycle=0; % the number of cycle can go below zero numerically. 
%     end;
%     
%     beta0=[beta0 n_components]; %Betti_0: the number of connected components
%     biggest0 = [biggest0 sizes(1)];  % the size of the largest component
%     beta1=[beta1 n_cycle];
%     
%     %temp(temp<rho)=0;   % remove edges below threshold rho
% end
% 
% biggest1=[]; %This is not yet implemnented. 


%---------------------------------------------------
% [nComponents,sizes,members] = networkComponents(A)
%
% Daniel Larremore
% April 24, 2014
% larremor@hsph.harvard.edu
% http://danlarremore.com
% Comments and suggestions always welcome.
% 
% INPUTS:
% A                     Matrix. This function takes as an input a 
% network adjacency matrix A, for a network that is undirected. If you
% provide a network that is directed, this code is going to make it
% undirected before continuing. Since link weights will not affect
% component sizes, weighted and unweighted networks work equally well. You
% may provide a "full" or a "sparse" matrix.
%
% OUTPUTS:
% nComponents             INT - The number of components in the network.
% sizes                 vector<INT> - a vector of component sizes, sorted, 
%   descending.
% members               cell<vector<INT>> a cell array of vectors, each
%   entry of which is a membership list for that component, sorted, 
%   descending by component size.
%
% Example: (uncomment and copy and paste into MATLAB command window)
% % Generate a 1000 node network adjacency matrix, A
% A = floor(1.0015*rand(1000,1000)); A=A+A'; A(A==2)=1; A(1:1001:end) = 0;
% % Call networkComponents function
% [nComponents,sizes,members] = networkComponents(A);
% % get the size of the largest component
% sizeLC = sizes(1);
% % get a network adjacency matrix for ONLY the largest component
% LC = A(members{1},members{1});

% function [nComponents,sizes,members] = networkComponents(A)
% % Number of nodes
% N = size(A,1);
% % Remove diagonals
% A(1:N+1:end) = 0;
% % make symmetric, just in case it isn't
% A=A+A';
% % Have we visited a particular node yet?
% isDiscovered = zeros(N,1);
% % Empty members cell
% members = {};
% % check every node
% for n=1:N
%     if ~isDiscovered(n)
%         % started a new group so add it to members
%         members{end+1} = n;
%         % account for discovering n
%         isDiscovered(n) = 1;
%         % set the ptr to 1
%         ptr = 1;
%         while (ptr <= length(members{end}))
%             % find neighbors
%             nbrs = find(A(:,members{end}(ptr)));
%             % here are the neighbors that are undiscovered
%             newNbrs = nbrs(isDiscovered(nbrs)==0);
%             % we can now mark them as discovered
%             isDiscovered(newNbrs) = 1;
%             % add them to member list
%             members{end}(end+1:end+length(newNbrs)) = newNbrs;
%             % increment ptr so we check the next member of this component
%             ptr = ptr+1;
%         end
%     end
% end
% % number of components
% nComponents = length(members);
% for n=1:nComponents
%     % compute sizes of components
%     sizes(n) = length(members{n});
% end
% 
% [sizes,idx] = sort(sizes,'descend');
% members = members(idx);

