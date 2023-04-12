function H = PH_hodge(B)
% PH_HODGE Compute all Hodge Laplacians of a simplicial complex.
%
%   H = PH_hodge(B) computes all Hodge Laplacians of a simplicial
%   complex given its boundary matrices B.
%
%   Inputs:
%       B - cell array containing the boundary matrices of the simplices.
%
%   Outputs:
%       H - cell array containing all Hodge Laplacians.
%
%
% Hodge Laplacian is first introduced in the medical imaging field in
%
% Lee, H., Chung, M.K., Kang, H., Lee, D.S. 2014. Hole detection in metabolic 
% connectivity of Alzheimer's disease using k-Laplacian, 17th International 
% Conference on Medical Image Computing and Computer Assisted Intervention (MICCAI). 
% Lecture Notes in Computer Science(LNCS) 8675:297-304
% https://pages.stat.wisc.edu/~mchung/papers/lee.2014.MICCAI.pdf
%
% The code is downloaded from
% https://github.com/laplcebeltrami/PH-STAT
%
% (C) 2023 Moo K. Chung, University of Wisconsin-Madison 
%     Email: mkchung@wisc.edu

% Get the number of boundary matrices
num_matrices = length(B);

% Initialize the cell array for Hodge Laplacians
H = cell(num_matrices, 1);

% Compute all Hodge Laplacians

k = 0;
H{k+1} = B{1} * B{1}'; % Graph Laplacian

for k = 1:(num_matrices - 1)
    H{k+1} =  B{k+1}*B{k+1}' + B{k}'*B{k};
end

