function betti = PH_hodge_betti(L)
%PH_HODGE_BETTI Compute the Betti numbers from Hodge Laplacians.
%
%   betti = PH_hodge_betti(L) computes the Betti numbers of a
%   simplicial complex given its Hodge Laplacians.
%
%   INPUT:
%       L - cell array of Hodge Laplacian matrices (can be sparse or full)
%
%   OUTPUT:
%       betti - vector of Betti numbers (betti(k) = dim ker L{k})
%
% (C) 2023 Moo K. Chung
%     University of Wisconsin–Madison
%
% Update: Works for sparse matrix as well, October 22, 2025


num_dim = length(L);
betti = zeros(1, num_dim);

for k = 1:num_dim
    H = full(L{k});                        % convert to full if sparse
    null_dim = size(H, 1) - rank(H);       % dim(ker(H)) = n - rank(H)
    betti(k) = null_dim;
end


% function betti = PH_hodge_betti(H_all)
% %   function betti = PH_hodge_betti(H_all) computes the Betti numbers of a
% %   simplicial complex given its Hodge Laplacians H_all.
% %
% %   Inputs:
% %       H_all - cell array containing the Hodge Laplacian of the simplices.
% %
% %   Outputs:
% %       betti - vector containing the Betti numbers.
% %
% % (C) 2023 Moo K. Chung, University of Wisconsin-Madison 
% %     Email: mkchung@wisc.edu
% 
% 
% % Get the number of Hodge Laplacians
% 
% %make it full matrix from sparse
% for i = 1:length(B)
%     B{i} = full(B{i});
% end
% 
% num_hodge = length(H_all);
% 
% % Initialize the Betti numbers vector
% betti = zeros(1, num_hodge);
% 
% % Compute the Betti numbers for each dimension
% for k = 1:num_hodge
%     H = H_all{k};
%     betti(k) = size(null(H),2);
% end