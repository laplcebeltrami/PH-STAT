function betti = PH_boundary_betti(B)
%PH_BOUNDARY_BETTI Compute the Betti numbers of a simplicial complex.
%
%   betti = PH_betti(B) computes the Betti numbers of the simplicial
%   complex represented by the boundary matrices in the cell array B.
%   The output is a vector containing the Betti numbers in each dimension.
%
%   Inputs:
%       B - cell array containing the boundary matrices of the simplices
%           in the complex.
%
%   Outputs:
%       betti - vector containing the Betti numbers of the complex. If
%       B{1}, B{2}, .. B{k} is available, up to betti{k-1} is given. 
%       betti{k} is the (k-1)-th Betti number. 
%
%
% (C) 2023 Moo K. Chung, Universtiy of Wisconsin-Madison 
%
%     Email: mkchung@wisc.edu
%
% Compute the number of dimensions of the simplicial complex
k = length(B);
p = size(B{1},1); %number of nodes

% Initialize the Betti numbers vector
betti = zeros(1,k);

%Betti-0
betti(1) = p - rank(B{1});

% Compute the Betti numbers for each dimension
for d = 2:k
    %(d-1)-th Betti number
    betti(d)=  rank(null(B{d-1})) - rank(B{d});
end

