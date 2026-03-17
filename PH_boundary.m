function B = PH_boundary(S)
% PH_BOUNDARY constructs boundary matrices from a simplicial complex.
%
% INPUT
%   S : cell array describing a simplicial complex.
%
%   Supported formats:
%
%   1) Pure simplex list
%      S{1} : vertices      (nV x 1)
%      S{2} : edges         (nE x 2)
%      S{3} : triangles     (nF x 3)
%      S{4} : tetrahedra    (nT x 4)
%
%   2) Two-column cell array
%      S{d,1} : simplex indices
%      S{d,2} : functional/scalar data attached to simplices
%
%      In this case, only S(:,1) is used here.
%
% OUTPUT
%   B : cell array of sparse integer boundary matrices
%
%       B{1} : boundary from 1-simplices to 0-simplices
%              rows = vertices, columns = edges
%
%       B{2} : boundary from 2-simplices to 1-simplices
%              rows = edges, columns = triangles
%
%       B{3} : boundary from 3-simplices to 2-simplices
%              rows = triangles, columns = tetrahedra
%
% IMPORTANT
%   The boundary matrices follow the simplex indexing and orientation
%   exactly as given in S. No sorting or reordering is performed.
%   Thus, if an edge in S{2} is listed as [3 1], then that orientation
%   is used as-is. Likewise, higher-dimensional simplices inherit the
%   orientation given by the vertex ordering in S.
%
%   The alternating sign convention is
%
%      d[v1,...,v_{k+1}] = sum_{j=1}^{k+1} (-1)^(j-1) [v1,...,v_{j-1},v_{j+1},...,v_{k+1}].
%
% (C) 2023- Moo K. Chung
%     University of Wisconsin-Madison
%     mkchung@wisc.edu
%
% Update history
%   2023        created
%   2025 Aug 26 documentation updated
%   2026 Mar 16 simplified and accelerated  

% If S is given as (#levels x 2) cell array, 
% use the first column only.  
if size(S,2) == 2 && size(S,1) > 1
    S = S(:,1);  
end

nLevel = numel(S);

if nLevel < 2
    B = {};
    return;
end

B = cell(nLevel-1,1);

for d = 1:(nLevel-1)
    L = S{d};     % (d-1)-simplices: nLower x d
    U = S{d+1};   % d-simplices:     nUpper x (d+1)

    nLower = size(L,1);
    nUpper = size(U,1);

    if isempty(L) || isempty(U)
        B{d} = sparse(nLower, nUpper);
        continue;
    end

    if size(L,2) ~= d || size(U,2) ~= d+1
        B{d} = sparse(nLower, nUpper);
        continue;
    end

    % Build lookup from lower-dimensional simplex to row index.
    % Orientation/indexing follows S exactly.  <— fixed
    key2row = containers.Map('KeyType','char','ValueType','double');
    for r = 1:nLower
        key2row(simplex_key(L(r,:))) = r;
    end

    % Preallocate at most (d+1)*nUpper nonzero entries.  <— fixed
    I = zeros((d+1)*nUpper, 1);
    J = zeros((d+1)*nUpper, 1);
    V = zeros((d+1)*nUpper, 1);
    nz = 0;

    for c = 1:nUpper
        up = U(c,:);

        for j = 1:(d+1)
            face = up([1:j-1, j+1:end]);   % oriented face induced by U
            fkey = simplex_key(face);

            if isKey(key2row, fkey)
                nz = nz + 1;
                I(nz) = key2row(fkey);
                J(nz) = c;
                V(nz) = (-1)^(j-1);
            end
        end
    end

    B{d} = sparse(I(1:nz), J(1:nz), V(1:nz), nLower, nUpper);
end

end

function k = simplex_key(v)
% SIMPLEX_KEY returns a row key without changing vertex order.
% Example: [3 1 5] -> '3,1,5'
k = sprintf('%d,', v);
k(end) = [];
end
