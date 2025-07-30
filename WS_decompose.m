function [Wb, Wd] = WS_decompose(W)
%WS_DECOMPOSE  Birth–death edge decomposition of weighted graphs.
%
%   [Wb, Wd] = WS_DECOMPOSE(W) returns the birth and death edge sets 
%   from a stack of weighted undirected adjacency matrices.
%
%   INPUT
%      W : edge weight matrix
%           - Can be either:
%           (i) a single p×p symmetric adjacency matrix (for one graph), or
%          (ii) a 3D array of size p×p×n, where n is the number of graphs (subjects).
%           - W(i,j,k) represents the weight of the edge between node i and j in the k-th graph.
%           - W is assumed to be symmetric with zero diagonal (no self-loops).
%
%   OUTPUT
%     Wb : (p–1) × 3 × n array of birth edges
%          Each slice Wb(:,:,k) contains [i, j, w] triplets of the first
%          (p–1) edges (by weight) that form the maximum spanning tree for subject k.
%
%     Wd : (p–1)(p–2)/2 × 3 × n array of death edges
%          Each slice Wd(:,:,k) contains the remaining edges not in the 
%          maximum spanning tree for subject k, ordered by descending weight.
%
%   EXAMPLE
%     For subject 1:
%        Wb(:,:,1) might be:
%           14   19   0.6314
%           20   22   0.6334
%            1    2   0.7321
%           ...
%        Wd(:,:,1) might be:
%            6   10   0.9375
%           10   11   0.8641
%            8   11   0.8513
%           ...
%
%   NOTE
%     - The function assumes input matrices are symmetric with zero
%     diagonals for fully connected complete graph
%     - Edges are sorted by weight, and MST (birth) edges are determined via Kruskal's algorithm.
%
%
% The method is published in
% 
% [1] Songdechakraiwut, T., Shen, L., Chung, M.K. 2021 Topological learning and 
% its application to multimodal brain network integration, Medical Image 
% Computing and Computer Assisted Intervention (MICCAI), LNCS 12902:166-176 
%
% [2] Songdechakraiwut, T. Chung, M.K. 2020 Topological learning for brain 
% networks, arXiv: 2012.00675. 
% 
% [3] Anand, D.V., Dakurah, S., Wang, B., Chung, M.K. 2021 Hodge-Laplacian 
% of brain networks and its application to modeling cycles. arXiv:2110.14599
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
%   2020 May 23 created. Modified from codes written by Song and Lee
%   2021 Nov 28 additional documentation
%   2022 Nov 22 conncomp_birth function is added
%   2024 Jul 25 For size(W,3) == 1, 
%               size is limited to 1000
%   2025 Mar 6  case size(birthMtx1, 1) < 100 and
%               size(deathMtx1, 1) < 2000 fixed
%   2026 Jul 28 Removed the size restriction, Added more comments for
%               PH-STAT manual. Loop simplified. 
%

[p, ~, nGraphs] = size(W);

% Preallocate outputs
Wb = zeros(p-1, 3, nGraphs);
Wd = zeros((p-1)*(p-2)/2, 3, nGraphs);

for k = 1:nGraphs
    Wk = W(:,:,k);

    % Compute MST (birth edges)
    birthEdges = conncomp_birth(Wk);
    Wb(:,:,k) = birthEdges;

    % Compute death edges
    G = graph(Wk, 'upper', 'omitselfloops');
    deathEdges = rmedge(G, birthEdges(:,1), birthEdges(:,2)).Edges{:, :};
    deathEdges = sortrows(deathEdges, 3, 'descend');
    Wd(:,:,k) = deathEdges;
end

% 
% if size(W,3) == 1
% 
%     % Compute set of births and set of deaths
%     G1 = graph(W, 'upper', 'omitselfloops');
%     %
%     % % birth edge set
%     birthMtx1 = conncomp_birth(W);
% 
%     % if size(birthMtx1, 1) > 100
%     %     indices = round(linspace(1, size(birthMtx1, 1), 100));
%     %     birthMtx2 = birthMtx1(indices, :);
%     % else
%     %     birthMtx2 = birthMtx1;
%     % end
% 
%     Wb = birthMtx1;
% 
%     %
%     % % death edge set
%     deathMtx1 = rmedge(G1, birthMtx1(:, 1), birthMtx1(:, 2)).Edges{:, :};
%     % % sorting by weights in ascending order
%     deathMtx1 = sortrows(deathMtx1, 3, 'descend');
% 
%     % if size(deathMtx1, 1) > 2000
%     %     indices = round(linspace(1, size(deathMtx1, 1), 1000));
%     %     deathMtx2 = deathMtx1(indices, :);
%     % else
%     %     deathMtx2 = deathMtx1;
%     % end
%     Wd = deathMtx1;
% 
% else %if there are more than 1 graph
%     %Wb = zeros(p-1, 3 , n);
%     % zeros((p-1)*(p-2)/2, 3, n);
% 
%     for i=1:size(W,3)
% 
%         Wi= W(:,:,i);
%         % Compute set of births and set of deaths
%         G1 = graph(Wi, 'upper', 'omitselfloops');
% 
%         % birth edge set
%         birthMtx1 = conncomp_birth(Wi);
%         Wb(:,:,i)=birthMtx1;
% 
%         % death edge set
%         deathMtx1 = rmedge(G1, birthMtx1(:, 1), birthMtx1(:, 2)).Edges{:, :};
%         % sorting by weights in ascending order
%         deathMtx1 = sortrows(deathMtx1, 3, 'descend');
% 
%         Wd(:,:,i)=deathMtx1;
% 
%     end
% end



function birthMtx = conncomp_birth(adj)
% Compute a set of increasing birth values for 0D barcode
%
% INPUT
%   adj      : weighted adjacency matrix
%
% OUTPUT
%   birthMtx : matrix whose 1st and 2nd columns are end nodes (no duplicates)
%              and 3rd column is weight (in ascending order, i.e., 1st row is
%              smallest)
%
% (C) 2020 Tananun Songdechakraiwut, Moo K. Chung
%          University of Wisconsin-Madison
%
%  Contact mkchung@wisc.edu for support with the codes 
%
% Update history:
%     2020 August 11 modified by Tananun from Lee's code
%     2021 May 25    comment by Chung

g = graph(-adj, 'upper', 'omitselfloops'); % minus weights to find max spanning tree
gTree = minspantree(g); % find max spanning tree of -adj
gTreeMtx = gTree.Edges{:, :}; % edge info.
gTreeMtx(:, 3) = gTreeMtx(:, 3) * -1; % reverse back to positive weights
birthMtx = sortrows(gTreeMtx, 3);



