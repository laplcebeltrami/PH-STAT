function cellcomplex_display_flows_colored(S, flow, vmin0, vmax0)
% CELLCOMPLEX_DISPLAY_FLOWS_COLORED Display directed edge flows with |flow| color.
% Optional:
%   vmin0, vmax0 : fixed color limits (override auto scaling)
%
% Usage:
%   cellcomplex_display_flows_colored(S, flow)
%   cellcomplex_display_flows_colored(S, flow, vmin0, vmax0)


% -------------------- parse Nodes --------------------
if isstruct(S) && isfield(S,'nodes')
    Nodes = S.nodes;
elseif iscell(S)
    Nodes = S{1};
else
    Nodes = S;
end

% accept flow vector aligned to S.edges
if ~isempty(flow) && size(flow,2)==1
    if isstruct(S) && isfield(S,'edges') && size(S.edges,1)==size(flow,1)
        flow = [S.edges, flow];
    else
        return
    end
end

if isempty(flow) || size(flow,2)~=3
    return
end

% basic sanity
K = size(Nodes,1);
if any(flow(:,1)<1) || any(flow(:,2)<1) || any(flow(:,1)>K) || any(flow(:,2)>K)
    return
end

hold on;
plot(Nodes(:,1), Nodes(:,2), 'k.', 'MarkerSize', 14);
axis equal;

% geometric parameters
delta   = 0.10;
shrink  = 0.22;
headL   = 0.18;
headW   = 0.1;

thr = 1e-8;

% -------------------- colormap scaling --------------------
magAll = abs(flow(:,3));
magAll = magAll(magAll > thr);

if isempty(magAll)
    hold off;
    return
end

% apply fixed limits if provided
if nargin >= 3 && ~isempty(vmin0) && ~isempty(vmax0)
    vmin = vmin0;
    vmax = vmax0;
else
    vmin = min(magAll);
    vmax = max(magAll);
end

if vmax == vmin
    vmax = vmin + eps;
end

% color setup
OrRd = brewermap(1000,'YlOrRd');
cmap = colormap(gca, OrRd(100:end,:));
colormap(gca, cmap);

caxis([vmin vmax]);
cb = colorbar;
cb.FontSize = 18;
cb.LineWidth = 1;

pos = cb.Position;
pos(2) = pos(2) + 0.5*(1 - 0.5)*pos(4);
pos(4) = 0.5 * pos(4);
cb.Position = pos;

% -------------------- draw edges --------------------
for k = 1:size(flow,1)

    i   = flow(k,1);
    j   = flow(k,2);
    val = flow(k,3);

    if abs(val) <= thr || i == j
        continue
    end

    if val > 0
        tail = i; head = j;
        mag  = val;
    else
        tail = j; head = i;
        mag  = -val;
    end

    u = (mag - vmin)/(vmax-vmin);
    u = max(0,min(1,u));
    idx = 1 + floor(u * (size(cmap,1)-1));
    col = cmap(idx,:);

    p1 = Nodes(tail,:);
    p2 = Nodes(head,:);

    d = p2 - p1;
    L = norm(d);
    if L == 0
        continue
    end

    tdir = d/L;
    ndir = [-tdir(2) tdir(1)];

    % shift & shrink
    p1 = p1 + delta*ndir + shrink*tdir;
    p2 = p2 + delta*ndir - shrink*tdir;

    pS = p2 - headL*tdir;

    lw = 3 + 1.8*(mag/(mag+1));

    plot([p1(1) pS(1)], [p1(2) pS(2)], 'Color', col, 'LineWidth', lw);

    v1 = p2;
    v2 = pS + headW*ndir;
    v3 = pS - headW*ndir;

    patch([v1(1) v2(1) v3(1)], ...
          [v1(2) v2(2) v3(2)], ...
          col, 'EdgeColor','none');
end

axis tight;
box off; axis off;
hold off;

end


% function cellcomplex_display_flows_colored(S, flow)
% % CELLCOMPLEX_DISPLAY_FLOWS_COLORED  Display directed edge flows with |flow| color.
% %
% % INPUT
% %   S    : struct with S.nodes (K x 2) and optionally S.edges (E x 2),
% %          or cell with S{1}=Nodes, or Nodes directly.
% %   flow : either
% %          (A) [E x 3] = [i j val]  (directed edge i->j with signed val), or
% %          (B) [E x 1] signed values aligned with S.edges (same E)
% %
% % COLOR / DIRECTION
% %   Color encodes |val|. Arrow direction encodes sign via orientation:
% %     if val>0, draw i -> j
% %     if val<0, draw j -> i
% %
% % WHEN EDGES ARE NOT DISPLAYED (common reasons)
% %   (1) abs(val) <= thr  (thresholded out)
% %   (2) val == 0         (no arrow drawn)
% %   (3) i==j             (self-loop ignored)
% %   (4) i or j out of range (not in 1..K)
% %   (5) Nodes are 2D but you pass 3D coordinates or vice versa (layout mismatch)
% %   (6) In case (B), length(flow) must match size(S.edges,1)
% %
% % (C) 2026 Moo K. Chung
% 
% % -------------------- parse Nodes --------------------
% if isstruct(S) && isfield(S,'nodes')
%     Nodes = S.nodes;
% elseif iscell(S)
%     Nodes = S{1};
% else
%     Nodes = S;
% end
% 
% % -------------------- accept flow as vector aligned with S.edges --------------------
% if ~isempty(flow) && size(flow,2) == 1
%     if isstruct(S) && isfield(S,'edges') && size(S.edges,1) == size(flow,1)
%         flow = [S.edges, flow];                    % [i j val]
%     else
%         return                                     % cannot infer edge indices
%     end
% end
% 
% if isempty(flow) || size(flow,2) ~= 3
%     return
% end
% 
% % -------------------- basic sanity (silent) --------------------
% K = size(Nodes,1);
% if any(flow(:,1)<1) || any(flow(:,2)<1) || any(flow(:,1)>K) || any(flow(:,2)>K)
%     return
% end
% 
% hold on;
% plot(Nodes(:,1), Nodes(:,2), 'k.', 'MarkerSize', 14);
% axis equal;
% 
% % -------------------- geometric parameters --------------------
% delta   = 0.10;
% shrink  = 0.22;
% headL   = 0.18;
% %headW   = 0.05;
% headW   = 0.1;
% 
% thr = 1e-8;  % <— edges with abs(val) <= thr are not drawn
% 
% % -------------------- colormap scaling (use only edges that will be drawn) --------------------
% magAll = abs(flow(:,3));
% magAll = magAll(magAll > thr);
% 
% if isempty(magAll)
%     hold off;
%     return
% end
% 
% vmin = min(magAll);
% vmax = max(magAll);
% if vmax == vmin
%     vmax = vmin + eps;
% end
% 
% %OrRd=brewermap(1000,'OrRd');
% OrRd=brewermap(1000,'YlOrRd');
% cmap = colormap(gca, OrRd(100:end,:));
% %cmap = colormap(gca, flipud(hot(256)));
% colormap(gca, cmap);  
% caxis([vmin vmax]);
% cb = colorbar;
% cb.FontSize = 18;
% cb.LineWidth = 1;
% 
% pos = cb.Position;
% pos(2) = pos(2) + 0.5*(1 - 0.5)*pos(4);
% pos(4) = 0.5 * pos(4);
% cb.Position = pos;
% 
% 
% % -------------------- draw edges --------------------
% for k = 1:size(flow,1)
% 
%     i   = flow(k,1);
%     j   = flow(k,2);
%     val = flow(k,3);
% 
%     % Not displayed if:
%     %   - abs(val) <= thr   (includes val==0 if thr==0)
%     %   - i==j              (self-loop)
%     if abs(val) <= thr || i==j
%         continue
%     end
% 
%     % direction from sign
%     if val > 0
%         tail = i; head = j;
%         mag  = val;
%     else
%         tail = j; head = i;
%         mag  = -val;
%     end
% 
%     % map magnitude to color
%     u   = (mag - vmin) / (vmax - vmin);
%     u   = max(0, min(1, u));
%     idx = 1 + floor(u * (size(cmap,1)-1));
%     col = cmap(idx,:);
% 
%     p1 = Nodes(tail,:);
%     p2 = Nodes(head,:);
% 
%     d = p2 - p1;
%     L = norm(d);
%     if L==0
%         continue
%     end
% 
%     tdir = d / L;
%     ndir = [-tdir(2) tdir(1)];
% 
%     % shift and shrink for cleaner arrows
%     p1 = p1 + delta*ndir;
%     p2 = p2 + delta*ndir;
% 
%     p1 = p1 + shrink*tdir;
%     p2 = p2 - shrink*tdir;
% 
%     % shaft endpoint (before arrowhead)
%     pS = p2 - headL*tdir;
% 
%     % linewidth (mild magnitude cue)
%     %lw = 1.2 + 1.8*(mag/(mag+1));
%     lw = 3 + 1.8*(mag/(mag+1));
% 
%     % shaft
%     plot([p1(1) pS(1)], [p1(2) pS(2)], 'Color', col, 'LineWidth', lw);
% 
%     % head
%     v1 = p2;
%     v2 = pS + headW*ndir;
%     v3 = pS - headW*ndir;
% 
%     patch([v1(1) v2(1) v3(1)], ...
%           [v1(2) v2(2) v3(2)], ...
%           col, 'EdgeColor','none');
% end
% 
% %figure_bg('w')          % requires your helper; if missing, comment out
% axis tight;
% box off; axis off;
% hold off;
% 
% end
