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

% -----------
% helper
% function draw_arrow_v(p1, p2, w, wmin, wmax, cmap)
% % draw_arrow_v  Draw a V-shaped arrow between two points with color mapping
% %
% % INPUT
% %   p1, p2 : 1x3 vectors (start and end node coordinates)
% %   w      : edge weight
% %   wmin   : minimum weight (for color normalization)
% %   wmax   : maximum weight (for color normalization)
% %   cmap   : colormap (Nx3)
% %
% % (C) 2026 Moo K. Chung
% % University of Wisconsin-Madison
% 
% d = p2 - p1;
% L = norm(d);
% if L == 0
%     return
% end
% 
% tdir = d / L;
% ndir = [-tdir(2), tdir(1), 0];
% nd = norm(ndir);
% if nd > 0
%     ndir = ndir / nd;
% else
%     ndir = [0 1 0];
% end
% 
% % Geometry
% headL = 0.2;
% headW = 0.1;
% 
% % Arrow placement (middle of edge)
% pm  = 0.5 * (p1 + p2);
% tip = pm + 0.18 * tdir;
% pS  = tip - headL * tdir;
% 
% % Line width
% if wmax == 0
%     wmax = 1;  % avoid division by zero
% end
% lw = 2 + 2 * w / wmax;
% 
% % Map weight to color
% if wmax == wmin
%     ucol = 0.5;
% else
%     ucol = (w - wmin) / (wmax - wmin);
% end
% ucol = max(0, min(1, ucol));
% idx = 1 + floor(ucol * (size(cmap,1)-1));
% col = cmap(idx,:);
% 
% % V-shaped arrow head
% v1 = tip;
% v2 = pS + headW * ndir;
% v3 = pS - headW * ndir;
% 
% plot3([v1(1) v2(1)], ...
%       [v1(2) v2(2)], ...
%       [v1(3) v2(3)], ...
%       'Color', col, ...
%       'LineWidth', lw);
% 
% plot3([v1(1) v3(1)], ...
%       [v1(2) v3(2)], ...
%       [v1(3) v3(3)], ...
%       'Color', col, ...
%       'LineWidth', lw);
% 
%   % Edge weight label near arrow head
%     pText = tip + 0.2 * ndir;
%     text(pText(1), pText(2), pText(3), ...
%          sprintf('%.2f', w), ...
%          'FontSize', 12, ...
%          'Color', col, ...
%          'HorizontalAlignment', 'center', ...
%          'VerticalAlignment', 'middle');
% 
% end