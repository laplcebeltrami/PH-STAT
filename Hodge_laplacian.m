function HLmat = Hodge_laplacian(IncidenceMat)
% function HLmat = Hodge_laplacian(IncidenceMat)
%
% Compute the 1st Hodge Laplacian for a graph using the incidence matrix.
%
% INPUT
%   IncidenceMat : cell array of boundary matrices from PH_boundary
%
% OUTPUT
%   HLmat        : 1st Hodge Laplacian matrix
%
% (C) 2021 Vijay Anand, Moo K. Chung
%     University of Wisconsin-Madison
%
% Update history
%   2021 November 11 created
%   2021 November 25 commented
%   2026 March 16 fixed indexing  

IndMat = IncidenceMat{1};     
HLmat  = IndMat' * IndMat;
end
