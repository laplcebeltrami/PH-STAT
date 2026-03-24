function Yvec = Hodge_vec(adj)
% Hodge_vec vectorizes the upper triangle of an adjacency matrix,
% excluding zero entries, in lexicographic order (i<j).
%
% INPUT
%   adj : n x n adjacency / connectivity matrix
%
% OUTPUT
%   Yvec : vector of length n*(n-1)/2 containing
%          adj(1,2), adj(1,3), ..., adj(1,n),
%          adj(2,3), ..., adj(n-1,n)
%
% (C) 2024 Vijay Anand, Moo K. Chung
% University of Wisconsin-Madison
%
% Update history; 2026 Mar 22 simplified

n = size(adj,1);

edges = find(triu(abs(adj) > 0, 1));
m = length(edges);
Yvec = zeros(m,1);

for e = 1:m
    [i,j] = ind2sub([n,n], edges(e));
    Yvec(e) = adj(i,j);
end

end