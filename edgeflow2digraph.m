function Gdir = edgeflow2digraph(kSkeleton, X)
% edgeflow2digraph  Convert edge flow into MATLAB digraph using kSkeleton.
%
% INPUT
%   kSkeleton : cell structure from Hodge_2Skeleton
%   X         : m x 1 edge flow vector aligned with kSkeleton{1,2}
%
% OUTPUT
%   Gdir      : MATLAB digraph object
%
% (C) 2026 Moo K. Chung
% University of Wisconsin-Madison

X = X(:);

elist = kSkeleton{1,2};        % <--- FIX: extract edge list
n = max(elist(:));             % <--- FIX: number of nodes
m = size(elist,1);

u = zeros(m,1);
v = zeros(m,1);
w = zeros(m,1);

for e = 1:m
    i = elist(e,1);
    j = elist(e,2);

    if X(e) >= 0
        u(e) = i;
        v(e) = j;
        w(e) = X(e);
    else
        u(e) = j;
        v(e) = i;
        w(e) = -X(e);
    end
end

idx = w > 1e-12;   % <--- FIX: remove zero edges

Gdir = digraph(u(idx), v(idx), w(idx), n);

end