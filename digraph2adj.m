function A = digraph2adj(G)
% Convert digraph to skew-symmetric adjacency matrix
%
% (C) 2026 Moo K. Chung
% University of Wisconsin-Madison

n = numnodes(G);
E = G.Edges.EndNodes;

if ismember('Weight', G.Edges.Properties.VariableNames)
    W = G.Edges.Weight(:);
else
    W = ones(size(E,1),1);
end

A = zeros(n,n);

for k = 1:size(E,1)
    i = E(k,1);
    j = E(k,2);

    A(i,j) = W(k);
    A(j,i) = -W(k);   
end

end