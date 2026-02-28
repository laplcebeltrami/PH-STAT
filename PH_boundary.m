function B = PH_boundary(S)
% PH_boundary constructs boundary matrices B from the augmented simplicial complex S
% where each S{d} has the last column as functional data.
% indexing in S follows standard orientation i<j, i<j<k
%
% INPUT:
%   S - Cell array: each S{d} is a [n x (d+1)] matrix with last column as
%   scalar data over simplices
%
%   S{1} : vertices       (nV x 1)          [optional; not used by many pipelines]
%   S{2} : edges          (nE x 2)  with i<j
%   S{3} : triangles      (nF x 3)  with i<j<k
%   S{4} : tetrahedra     (nT x 4)  with i<j<k<l
%
% OUTPUT:
%  B - Cell array of sparse integer boundary matrices B{d}, 
%       each mapping (d-1)-simplices to d-simplices
%
%       B{1}: The vertex-to-edge boundary matrix.
%           - Rows: vertices (0-simplices)
%           - Columns: edges (1-simplices)
%           - An entry B{1}(i,j) is nonzero if vertex i is a face of edge j.
%             The sign (+1 or -1) is determined by the ordering (orientation)
%             of the edge.
%
%       B{2}: The edge-to-triangle boundary matrix.
%           - Rows: edges (1-simplices)
%           - Columns: triangles (2-simplices)
%           - An entry B{2}(i,j) is nonzero if edge i is a face of triangle j.
%             The sign is determined by the orientation of the triangle.
%
%       B{3}: (If available) The triangle-to-tetrahedron boundary matrix.
%           - Rows: triangles (2-simplices)
%           - Columns: tetrahedra (3-simplices)
%           - An entry B{3}(i,j) is nonzero if triangle i is a face of 
%             tetrahedron j, with the sign reflecting the orientation.
%
%   Note 1. The function uses the alternating sign convention:
%   For a d-simplex [v1, v2, ..., v(d+1)], its boundary is given by
%   ∂[v1, v2, ..., v(d+1)] 
%   = Σ_{j=1}^{d+1} (-1)^(j-1) [v1,...,v(j-1),v(j+1),...,v(d+1)].
%
%   Note 2. B{1}' assumes columns of B{1} are edges in reference orientation i<j.
%           B{2}' assumes columns of B{2} are triangles in reference orientation i<j<k.
%
% (C) 2023- Moo K. Chung
%     University of Wisconsin-Madison
%     Email: mkchung@wisc.edu
%
%     Update histroy: 2023 created
%                    2025 August 26 updated with documents


 nC = numel(S);
    if nC < 2
        B = {};
        return;
    end

    B = cell(nC-1,1);

    for d = 1:(nC-1)
        % pull levels (may be empty or too narrow)
        L = S{d};       % expect n x d
        U = S{d+1};     % expect n x (d+1)

        % Determine sizes safely
        nLower = 0; nUpper = 0;
        if ~isempty(L), nLower = size(L,1); end
        if ~isempty(U), nUpper = size(U,1); end

        % If either level is missing, return correctly-sized sparse zeros
        if isempty(L) || isempty(U) || size(L,2) < d || size(U,2) < (d+1)
            B{d} = sparse(nLower, nUpper);
            continue;
        end

        % Keep only the vertex columns and sort to canonical orientation
        L = sort(L(:,1:d),   2);   % faces (rows) in reference orientation
        U = sort(U(:,1:d+1), 2);   % simplices (cols) in reference orientation

        % Build face->row map
        face2row = containers.Map('KeyType','char','ValueType','double');
        for r = 1:size(L,1)
            face2row( key_of(L(r,:)) ) = r;
        end

        % Fill boundary matrix using alternating sign convention
        Bmat = sparse(size(L,1), size(U,1));
        for c = 1:size(U,1)
            up = U(c,:);
            for idx = 1:(d+1)
                face = up([1:idx-1, idx+1:end]);  % drop vertex idx
                fkey = key_of(face);
                if isKey(face2row, fkey)
                    r = face2row(fkey);
                    Bmat(r,c) = (-1)^(idx-1);
                end
            end
        end

        B{d} = Bmat;
    end
end

% ---- helper: canonical key for a row of indices ----
function k = key_of(v)
    % e.g., [3 7 10] -> '3,7,10'
    k = sprintf('%d,', v);
    k(end) = [];
end

% 
% k = length(S) - 1;         % Max simplex dimension
% B = cell(k, 1);            % Allocate output
% 
% for d = 1:k
%     lowerSimplices = S{d}(:,1:d);      % Drop scalar column, keep (d)-vertex list
%     upperSimplices = S{d+1}(:,1:d+1);  % Drop scalar column, keep (d+1)-vertex list
% 
%     nLower = size(lowerSimplices, 1);
%     nUpper = size(upperSimplices, 1);
%     Bmat = sparse(nLower, nUpper);
% 
%     % Build map from face to row index
%     keyMap = containers.Map;
%     for i = 1:nLower
%         face = sort(lowerSimplices(i,:));
%         key = join(string(face), ',');
%         keyMap(key) = i;
%     end
% 
%     % Fill in boundary matrix
%     for j = 1:nUpper
%         simplex = upperSimplices(j,:);
%         for idx = 1:(d+1)
%             face = simplex;
%             face(idx) = [];                  % Remove one vertex
%             face = sort(face);               % Sort to standardize
%             key = join(string(face), ',');
%             if isKey(keyMap, key)
%                 row = keyMap(key);
%                 Bmat(row, j) = (-1)^(idx - 1);
%             end
%         end
%     end
% 
%     B{d} = Bmat;
% end
% end
% 
% 
% 
% 
