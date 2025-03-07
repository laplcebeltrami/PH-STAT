function S = PH_connectivity2simplex(C, tau, maxDim)
% PH_connectivity2simplex constructs a simplicial complex from a 
% correlation matrix using a single threshold tau, building higher-dimensional
% simplices iteratively from lower-dimensional ones.
%
%   S = PH_connectivity2simplex(C, tau, maxDim) builds a simplicial
%   complex from the weighted connectivity matrix C (p x p) as follows:
%
%       - Nodes (0-simplices) are the p time series (or brain regions).
%
%       - An edge (1-simplex) between nodes i and j is included if C(i,j) > tau.
%
%       - A triangle (2-simplex) among nodes i, j, and k is included if every pair
%         among {i, j, k} has C > tau.
%
%       - A tetrahedron (3-simplex) among nodes i, j, k, and l is included if every
%         pair among {i, j, k, l} has C > tau.
%
%       - And so on, up to the maximum simplex dimension maxDim.
%
%   The function returns a cell array S, where:
%       S{1} is a column vector of nodes,
%       S{2} is a matrix of edges (each row is a 2-element vector),
%       S{3} is a matrix of triangles (each row is a 3-element vector),
%       S{4} is a matrix of tetrahedra (each row is a 4-element vector), etc.
%
%   Inputs:
%       C      - p x p correlation matrix (assumed symmetric, with ones on the diagonal)
%       tau    - scalar threshold used uniformly to include simplices
%       maxDim - maximum dimension of the simplicial complex 
%                (e.g., maxDim = 3 builds up to tetrahedra, so S will have 4 cells)
%
%   Output:
%       S      - Cell array representing the simplicial complex.
%
% (C) 2025 Moo K. Chung
%   University of Wisconsin-Madison
%     Email: mkchung@wisc.edu

% Get the number of nodes.
p = size(C,1);

%% Build the Simplicial Complex S Iteratively

% S{1}: 0-simplices (nodes)
% Each node is simply represented by its index.
S{1} = (1:p)';

% S{2}: 1-simplices (edges)
% Include each edge (i,j) with i < j if C(i,j) > tau.
edges = [];
for i = 1:p-1
    for j = i+1:p
        if C(i,j) > tau
            edges = [edges; i, j]; %#ok<AGROW>
        end
    end
end
S{2} = edges;

% For higher dimensions d>=2, build S{d+1} from S{d} iteratively.
for d = 2:maxDim
    % S{d} contains the (d-1)-simplices.
    prevSimplices = S{d};
    newSimplices = [];
    
    % Iterate over each (d-1)-simplex.
    for i = 1:size(prevSimplices, 1)
        simplex = prevSimplices(i, :);
        % To avoid duplicates, we only add a new vertex v that is greater than the maximum
        % vertex in the current simplex.
        m = max(simplex);
        for v = m+1:p
            % Check if vertex v is connected to every vertex in 'simplex'
            % i.e., for every u in 'simplex', we need C(v,u) > tau.
            if all(C(v, simplex) > tau)
                % If so, the new simplex is the union of the current simplex and vertex v.
                newSimplex = sort([simplex, v]);  % Sorting is optional if simplex is already sorted.
                newSimplices = [newSimplices; newSimplex]; %#ok<AGROW>
            end
        end
    end
    
    % Store the new simplices in S{d+1}. 
    % If no new simplices are found, S{d+1} will be an empty matrix.
    S{d+1} = newSimplices;
end

end