function S = PH_graphfiltration(X, e)
%   S = PH_graphfiltration(X,e) computes the 1-skeleton out of graph
%   filtration at filtration value e
%
%   Inputs:
%       X - a p x d matrix of coordinates of scatter points in R^d
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
maxw = max(w(:)); %maximum edge weight
e=maxw-e; %Graph filtraiton is equivalent to Rips filtration on 1-e
% Compute the number of points in the cloud
p = size(w, 1);

% Initialize the output cell array
k=1
S = cell(k+1, 1);

% Add vertice coordinates to the output cell array
S{1} = (1:p)';


% Add higher-dimensional simplices to the output cell array

d=2
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


