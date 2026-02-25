function vec = adj2vec(adj, k);
%
% function vec = adj2vec(adj,k);
% vectorize the square matrix and produce the vector of elemenets 
% in the upper triangle. For example, from the distance matrix, it produces the pdist.
% 
% INPUT:
% adj has to be bigger than 2 x 2 matrix
% k   row index. If k=1. It vectorize from the diagonal. 
% If k=2, it vectorize above the 2nd row exluding diagonal.
%
% Given
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
% mkchung@wics.edu
%


n=size(adj,1);
vec=adj(1,2:end);
 
if n>=2
    for i=k:n
        vec= [vec adj(i,i+1:end)];
    end; 
end
