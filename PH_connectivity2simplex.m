function S = PH_connectivity2simplex(C, tau, maxDim)
% PH_connectivity2simplex constructs a simplicial complex from a
% correlation matrix using threshold tau. Tries to speed up the building 
% of higher-dimensional simplices via adjacency lists and set intersections.
%
%   S = PH_connectivity2simplex(C, tau, maxDim) builds a simplicial
%   complex from the weighted connectivity matrix C (p x p) as follows:
%
%       - Nodes (0-simplices) are the p time series (or brain regions).
%       - An edge (1-simplex) between nodes i and j is included if C(i,j) > tau.
%       - A triangle (2-simplex) is included if each pair of vertices 
%         in that triple has C > tau, etc., up to maxDim.
%
%   The function returns a cell array S, where:
%       S{1} is a column vector of nodes,
%       S{2} is a matrix of edges,
%       S{3} is a matrix of triangles, etc.
%
%   Inputs:
%       C      - p x p correlation matrix (assumed symmetric, ones on diagonal)
%       tau    - threshold for including simplices
%       maxDim - maximum dimension of the simplicial complex 
%
%   Output:
%       S      - Cell array representing the simplicial complex.
%
%  The function is downloaded from https://github.com/laplcebeltrami/PH-STAT
%  and part of PH-STAT pacakge.
%
% (C) 2025 Moo K. Chung
%     University of Wisconsin-Madison
%     Email: mkchung@wisc.edu

p = size(C,1);
S{1} = (1:p)';

% Build adjacency matrix.
A = (C > tau);

% Build edges (1-simplices) using the upper triangle of A.
[i,j] = find(triu(A,1));
S{2} = [i,j];

% Create adjacency lists: adjacencyList{n} = row vector of neighbors of node n.
% Using cell arrays for flexible storage.
adjacencyList = cell(p,1);
for n = 1:p
    adjacencyList{n} = find(A(n,:));
end

% Build higher-dimensional simplices up to maxDim.
for d = 2:maxDim
    prevSimplices = S{d};
    newSimplices = [];
    
    for k = 1:size(prevSimplices,1)
        simplex = prevSimplices(k,:);
        m = max(simplex);
        
        % Intersect adjacency lists of each node in 'simplex' to find common neighbors.
        % Start with the neighbors of the first node, then intersect repeatedly.
        commonNeighbors = adjacencyList{simplex(1)};
        for idx = 2:length(simplex)
            commonNeighbors = intersect(commonNeighbors, adjacencyList{simplex(idx)});
            if isempty(commonNeighbors)
                break;  % Faster exit if no intersection remains
            end
        end
        
        % Keep only those neighbors that exceed 'm' to avoid duplicates.
        commonNeighbors = commonNeighbors(commonNeighbors > m);
        
        % Form new (d+1)-simplices.
        for v = commonNeighbors
            newSimplex = [simplex, v];
            newSimplices = [newSimplices; sort(newSimplex)]; %#ok<AGROW>
        end
    end
    
    % Remove duplicates (if any).
    if ~isempty(newSimplices)
        newSimplices = unique(newSimplices,'rows');
    end
    S{d+1} = newSimplices;
end

end
