function betti = PH_hodge_betti(H_all)
%   function betti = PH_hodge_betti(H_all) computes the Betti numbers of a
%   simplicial complex given its Hodge Laplacians H_all.
%
%   Inputs:
%       H_all - cell array containing the Hodge Laplacian of the simplices.
%
%   Outputs:
%       betti - vector containing the Betti numbers.
%
% (C) 2023 Moo K. Chung, University of Wisconsin-Madison 
%     Email: mkchung@wisc.edu


% Get the number of Hodge Laplacians
num_hodge = length(H_all);

% Initialize the Betti numbers vector
betti = zeros(1, num_hodge);

% Compute the Betti numbers for each dimension
for k = 1:num_hodge
    H = H_all{k};
    betti(k) = size(null(H),2);
end