function [varargout] = Hodge_decompose(Yvec, Bmat)
% Hodge_decompose computes the Hodge decomposition of a vectorized
% connectivity matrix (edge flow).
%
% INPUT
%   Yvec : edge flow vector (nEdge x 1) or  vectorized connectivity matrix (edge flow)
%          For a p x p connectivity matrix, length(Yvec) should be p*(p-1)/2.
%
%           The ordering of Yvec MUST exactly match the column ordering
%           of the boundary matrix Bmat{1}. Equivalently, Yvec(k) must
%           correspond to the k-th edge listed in the simplicial complex.
%
%           Example:
%
%               S{2} =
%                   1   2   1.0
%                   2   3   1.0
%                   1   3  -0.5
%
%           Then:
%
%               Edges =
%                   [1 2
%                    2 3
%                    1 3]
%
%               Yvec =
%                   [ 1.0
%                     1.0
%                    -0.5 ]
%
%           Here:
%
%               Yvec(1) corresponds to edge (1,2)
%               Yvec(2) corresponds to edge (2,3)
%               Yvec(3) corresponds to edge (1,3)
%
%           The sign encodes direction relative to the edge orientation.
%
%           If canonical orientation i<j is used:
%
%               Yvec(k) > 0  means flow i -> j
%               Yvec(k) < 0  means flow j -> i
%
%           Example:
%
%               edge = (1,4)
%               Yvec = -1
%
%           means directed flow:
%
%               4 -> 1
%
%   Bmat : boundary matrices from Hodge_incidence.m or PH_boundary.m
%
%           Bmat{1} : node-edge incidence matrix
%                      size = nNode x nEdge
%
%           Bmat{2} : edge-triangle incidence matrix
%                      size = nEdge x nTriangle
%
% OUTPUT
%   Yg : gradient flow
%   Yc : curl flow
%   Yh : harmonic flow
%   Optional:
%   s  : node potential
%   z  : face potential
%
% (C) 2022 Vijay Anand D, Moo K. Chung
% University of Wisconsin-Madison
%
% Update history: 2022 April 24 created; 2022 Nov 15 comment; 
% 2024 July 30 5-output option; 2024 Oct 4 sparse matrices; 
% 2026 Mar 21 fixed boundary indexing  
    
    Yvec = Yvec(:);   % force column vector

    % Bmat{1}: node-edge boundary matrix
    d0 = sparse(Bmat{1}');   % edge-node operator  

    % Bmat{2}: edge-triangle boundary matrix, if available
    if length(Bmat) >= 2
        d1 = sparse(Bmat{2}');   % triangle-edge operator  
    else
        d1 = [];
    end

    
    % Hodge Laplacian for gradient part
    L0 = d0' * d0;
    
    % Divergence
    div = -d0' * Yvec;

    % Solve for gradient component
    %s = lsqr(L0, -div); 
    [s, flag_s] = lsqr(L0, -div);   
    Yg = d0 * s;
    
    % Curl component
    if ~isempty(d1)
        L1 = d1 * d1';
        curl = d1 * Yvec;
        %z = lsqr(L1, curl);
           [z, flag_z] = lsqr(L1, curl); 
        Yc = d1' * z;
        Yc=full(Yc);
    else
        z = [];
        Yc = zeros(size(Yvec));   
    end

    % Harmonic component
    Yh = Yvec - Yg - Yc;

    % Outputs
    if nargout == 3
        varargout = {Yg, Yc, Yh};
    elseif nargout == 5
        varargout = {Yg, Yc, Yh, s, z};
    end
end

