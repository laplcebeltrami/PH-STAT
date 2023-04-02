function S = PH_rips(X, k, e)
%   S = PH_rips(X, k, e) computes the k-skeleton of the Rips complex
%   of a point cloud given scatter points X and a radius e. The output
%   is a cell array S containing the simplices of the k-skeletons.
%
%   Inputs:
%       X - a p x d matrix of coordinates of scatter points in R^d
%       k - the maximum dimension of simplices to include
%       e - the radius parameter for the Rips complex
%
%   Outputs:
%
%       S - structured array containing nodes, edges, faces, volumes, etc.
%
%
% (C) 2023 Moo K. Chung, Universtiy of Wisconsin-Madison 
%
%     Email: mkchung@wisc.edu


% Compute the pairwise distances between points
w = pdist2(X,X);

% Compute the number of points in the cloud
p = size(w, 1);

% Initialize the output cell array
S = cell(k+1, 1);

% Add vertice coordinates to the output cell array
S{1} = (1:p)';


% Add higher-dimensional simplices to the output cell array
for d = 2:(k+1)
    % We compute all possible d-simplices using nchoosek. The loop iterates over
    % all combinations of d vertices from the point cloud, and adds the simplex
    % to S if all pairwise distances between vertices are at most e.

    % Compute the list of possible d-simplices
    possible_simplices = nchoosek(1:p, d);
    num_possible_simplices = size(possible_simplices, 1);
    valid_simplices = [];
    for ii = 1:num_possible_simplices
        simplex = possible_simplices(ii,:);
        % Check if all pairwise distances are at most e
        pairwise_distances = w(simplex, simplex);
        if all(pairwise_distances(:) <= e)
            valid_simplices = [valid_simplices; simplex];
        end
    end
    % Add the valid d-simplices to S
    S{d} = valid_simplices;
end

end


