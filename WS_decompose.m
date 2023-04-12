function [Wb Wd] = WS_decompose(W)
%function [Wb Wd] = WS_decompose(W)
%    
% Performs the birth-death decomposition of edge sets.
%
% INPUT
% W : Connectivity matrix, edge weight matrix or weighted adjacency martrix
%     of size p x p
%
% OUTPUT
% Wb : birth edge set     (p-1) x 3 x n, where p is # of nodes and n is # of subjects
% Wd : death edge set     (p-1)*(p-2)/2 x 3 x n, where p is # of nodes and n is # of subjects
%
%
% The method is published in
% 
% [1] Songdechakraiwut, T., Shen, L., Chung, M.K. 2021 Topological learning and 
%its application to multimodal brain network integration, Medical Image 
%Computing and Computer Assisted Intervention (MICCAI), LNCS 12902:166-176 
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


if size(W,3) == 1

    % Compute set of births and set of deaths
    G1 = graph(W, 'upper', 'omitselfloops');
    %
    % % birth edge set
    birthMtx1 = conncomp_birth(W);
    Wb=birthMtx1;
    %
    % % death edge set
    deathMtx1 = rmedge(G1, birthMtx1(:, 1), birthMtx1(:, 2)).Edges{:, :};
    % % sorting by weights in ascending order
    deathMtx1 = sortrows(deathMtx1, 3, 'descend');
    Wd=deathMtx1;

else %if there are more than 1 graph
    %Wb = zeros(p-1, 3 , n);
    % zeros((p-1)*(p-2)/2, 3, n);
  
    for i=1:size(W,3)
        
        Wi= W(:,:,i);
        % Compute set of births and set of deaths
        G1 = graph(Wi, 'upper', 'omitselfloops');

        % birth edge set
        birthMtx1 = conncomp_birth(Wi);
        Wb(:,:,i)=birthMtx1;

        % death edge set
        deathMtx1 = rmedge(G1, birthMtx1(:, 1), birthMtx1(:, 2)).Edges{:, :};
        % sorting by weights in ascending order
        deathMtx1 = sortrows(deathMtx1, 3, 'descend');

        Wd(:,:,i)=deathMtx1;

    end
end





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



