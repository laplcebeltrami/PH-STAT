function vec = adj2vec(adj, k)
% function vec = adj2vec(adj, k)
%
% Vectorizes a square matrix by extracting upper triangular elements.
%
% INPUT
%   adj : square matrix (n x n)
%   k   : optional starting row index
%         k = 1 → include from first row (default)
%         k = 2 → skip first row (exclude diagonal behavior similar
%         to pdist)
%
% OUTPUT
%   vec : row vector of selected upper-triangular elements
%
% EXAMPLE
% adj = [ 0  1  2  3;
%        1  0  4  5;
%        2  4  0  6;
%        3  5  6  0 ];
%
% It outputs vec = [1 2 3 4 5 6]
%
%
% (C) 2019 Moo K. Chung 
% University of Wisconsin-Madison
%
% Update history
%   2026 Mar 17 unified optional argument handling 

if nargin < 2
    k = 1;   
end

n = size(adj,1);

vec = [];    

for i = k:n
    if i < n
        vec = [vec adj(i, i+1:end)];
    end
end

end


