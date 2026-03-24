function coord = digraph_plot(Gdir, coord)
% digraph_plot(Gdir)
%
% Plot a directed graph by first displaying the underlying undirected
% topology and then overlaying dominant directed arrows.
%
% INPUT
%   Gdir : MATLAB digraph object
%   coord : (optional) n x 2 or n x 3 node coordinates
%
% OUTPUT
%   none
%
% (C) 2026 Moo K. Chung
% University of Wisconsin-Madison

% Convert directed graph to adjacency matrix
A = digraph2adj(Gdir);
A_undir = max(abs(A), abs(A')); %do not use abs. 
% It doubles for undirected graph
% Build undirected graph from the absolute adjacency matrix
G = graph(A_undir);


%p = plot(G);        % automatic layout (force-directed by default)
%coord = [p.XData(:), p.YData(:)];   % extracted node coordinates

% Plot graph with or without provided coordinates
if nargin < 2 || isempty(coord)
    p = plot(G);                             % automatic layout
    coord = [p.XData(:), p.YData(:)];        % extracted coordinates
else
    if size(coord,2) == 2
        p = plot(G,'XData',coord(:,1),'YData',coord(:,2));  % use given 2D coords
    end
end

% Plot undirected topology (assumes plot_graph uses current axes)
graph_plot(G,coord);
hold on;

% Overlay dominant directed arrows
graph_add_arrows(Gdir,coord);
end


%---
% helper
function graph_add_arrows(Gdir, coord)
% graph_add_arrows(Gdir, coord)
%
% Add dominant directed arrows on top of an existing undirected plot_graph(G)
% figure. If coord is provided, node coordinates are taken directly from it.
% Otherwise, coordinates are recovered from the numeric node labels placed
% by plot_graph.
%
% INPUT
%   Gdir  : MATLAB digraph object
%   coord : optional n x 2 or n x 3 node coordinates
%
% (C) 2026 Moo K. Chung
% University of Wisconsin-Madison

ax = gca;
hold(ax,'on');

n = numnodes(Gdir);

% -------------------------------------------------------------------------
% Get node coordinates

X = coord(:,1);
Y = coord(:,2);

if size(coord,2) >= 3
    Z = coord(:,3);
else
    Z = zeros(n,1);   
end

% -------------------------------------------------------------------------
% Extract dominant directed flows

flows = digraph2edgeflow_dominant(Gdir);
elist = flows(:,1:2);
weights = flows(:,3);



% -------------------------------------------------------------------------
% Colormap setup for whole figure
cmap = hot(256);
cmap = cmap(end-70:-1:1, :);
cmap = interp1(linspace(0,1,size(cmap,1)), cmap, linspace(0,1,256));
colormap(gca, cmap);

wmin = min(weights);
wmax = max(weights);
if wmax == wmin
    wmax = wmin + eps;
end

caxis([wmin wmax]);
cb = colorbar;
cb.FontSize = 16;
cb.LineWidth = 1;
cb.TickLabels = compose('%.1f', cb.Ticks);

pos = cb.Position;
pos(2) = pos(2) + pos(4)*0.25;
pos(4) = pos(4)*0.5;
cb.Position = pos;

% Draw arrows only

for e = 1:size(elist,1)
    
    u = elist(e,1);
    v = elist(e,2);

    p1 = [X(u), Y(u), Z(u)];
    p2 = [X(v), Y(v), Z(v)];

    draw_arrow_v(p1, p2, weights(e), wmin, wmax, cmap);   % <--- fixed
end


end


function draw_arrow_v(p1, p2, w, wmin, wmax, cmap)
% draw_arrow_v  Draw a V-shaped arrow between two points with color mapping
%
% INPUT
%   p1, p2 : 1x3 vectors (start and end node coordinates)
%   w      : edge weight
%   wmin   : minimum weight (for color normalization)
%   wmax   : maximum weight (for color normalization)
%   cmap   : colormap (Nx3)
%
% (C) 2026 Moo K. Chung
% University of Wisconsin-Madison

d = p2 - p1;
L = norm(d);
if L == 0
    return
end

tdir = d / L;
ndir = [-tdir(2), tdir(1), 0];
nd = norm(ndir);
if nd > 0
    ndir = ndir / nd;
else
    ndir = [0 1 0];
end

% Geometry
headL = 0.2;
headW = 0.1;

% Arrow placement (middle of edge)
pm  = 0.5 * (p1 + p2);
tip = pm + 0.18 * tdir;
pS  = tip - headL * tdir;

% Line width
if wmax == 0
    wmax = 1;  % avoid division by zero
end
lw = 2 + 2 * w / wmax;

% Map weight to color
if wmax == wmin
    ucol = 0.5;
else
    ucol = (w - wmin) / (wmax - wmin);
end
ucol = max(0, min(1, ucol));
idx = 1 + floor(ucol * (size(cmap,1)-1));
col = cmap(idx,:);

% V-shaped arrow head
v1 = tip;
v2 = pS + headW * ndir;
v3 = pS - headW * ndir;

plot3([v1(1) v2(1)], ...
      [v1(2) v2(2)], ...
      [v1(3) v2(3)], ...
      'Color', col, ...
      'LineWidth', lw);

plot3([v1(1) v3(1)], ...
      [v1(2) v3(2)], ...
      [v1(3) v3(3)], ...
      'Color', col, ...
      'LineWidth', lw);

  % Edge weight label near arrow head
    pText = tip + 0.2 * ndir;
    text(pText(1), pText(2), pText(3), ...
         sprintf('%.2f', w), ...
         'FontSize', 12, ...
         'Color', col, ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle');

end
%-----
function flows = digraph2edgeflow_dominant(G)
% digraph2edgeflow_dominant  Extract dominant directed edge flows from a digraph.
%
% INPUT
%   G : digraph object in MATLAB format
%
% OUTPUT
%   flows : r x 3 array [u v w]
%           u -> v is the dominant direction
%           w >= 0 is the dominant flow magnitude
%
% For each unordered node pair {i,j}, the function combines the two
% directed edges i->j and j->i into a single dominant directed flow.
%
% (C) 2026 Moo K. Chung
% University of Wisconsin-Madison

A = full(adjacency(G, 'weighted'));
n = size(A,1);

flows = [];

for i = 1:n
    for j = i+1:n

        wij = A(i,j);
        wji = A(j,i);

        % net dominant flow
        w = wij - wji;

        if abs(w) < 1e-12
            continue
        end

        if w > 0
            u = i;
            v = j;
            mag = w;
        else
            u = j;
            v = i;
            mag = -w;
        end

        flows = [flows; u v mag]; %#ok<AGROW>
    end
end

end