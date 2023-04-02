function B = PH_boundary(S)
%PH_boundary compute the boundary matrices of a simplicial complex.
%
%   B = PH_boundary(S) computes the boundary matrices of the simplicial
%   complex S represented as a cell array of simplices. The output is a
%   cell array B containing the boundary matrices of the complex in each
%   dimension.
%
%   Inputs:
%       S - cell array containing the simplices of the k-skeletons of a
%           Rips complex. It contains list of nodes that forms the k-simplex.
%
%   Outputs:
%       B - cell array containing the boundary matrices of the simplices in S.
%
% (C) 2023 Moo K. Chung, Universtiy of Wisconsin-Madison 
%
%     Email: mkchung@wisc.edu


% Compute the number of dimensions of the simplicial complex
k = length(S) - 1;

% Initialize the boundary matrices
B = cell(k, 1);



p = size(S{1}, 1);

% Compute the number of simplices in each dimension
num_simplices = zeros(1, k+1);
for d = 0:k
    num_simplices(d+1) = nchoosek(p, d+1);
end

% Compute the boundary matrices for each dimension
for d = 1:k
    % Get the simplices of dimension d and d+1
    simplices_d = S{d};
    simplices_d1 = S{d+1};

    % Compute the boundary matrix for dimension d
    num_simplices_d = size(simplices_d,1);
    num_simplices_d1 = size(simplices_d1,1);
    %boundary = sparse(num_simplices_d, num_simplices_d1);
    boundary = zeros(num_simplices_d, num_simplices_d1);

    for ii = 1:num_simplices_d1
        simplex = simplices_d1(ii, :);
        for jj = 1:d+1
            face = simplex;
            face(jj) = [];
            idx = find(ismember(simplices_d, face, 'rows'));
            if ~isempty(idx)
                boundary(idx, ii) = (-1)^(jj-1);
            end
        end
    end
    B{d} = boundary;
end





