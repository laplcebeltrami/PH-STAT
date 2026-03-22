function h=cellcomplex_display(S)
% CELLCOMPLEX_DISPLAY  Display a 2D cell complex and color ALL cycles,
% including overlapping cycles that share edges.
%
% Key fix:
%   When two cycles share edges, the shared edge cannot be oriented
%   consistently for both cycles at once. If we first plot the graph and then
%   fill faces only once, one of the overlapping cycles may not appear as a
%   closed polygon under a single global edge orientation. The remedy is:
%   build the polygon for EACH cycle independently (using its own local
%   ordering) and patch it on top of existing patches. Overlaps are handled
%   naturally by layering with transparency.
%
% Input
%   S.nodes : [K x 2] node coordinates
%   S.edges : [E x 2] directed edge list
%
% Supported face encodings:
%   (A) S.faces.vertices : [F x q] vertex indices
%   (B) S.faces.cycles   : cell array; each entry is [L x 2] edge list
%
% (C) 2026 Moo K. Chung

Nodes = S.nodes;
Edges = S.edges;

% -------------------- base graph (draw once) --------------------
%edgecolor = [0.1 0.1 0.1]; % dark gray
edgecolor = [0.6 0.6 0.6];   % light gray
edgelw    = 20;

graph_directed_plot(Nodes, Edges, [], edgecolor);
hold on;

% enforce edge styling uniformly (graph_directed_plot might differ)
h = findobj(gca,'Type','Line');
set(h,'LineWidth',edgelw,'Color',edgecolor);

% -------------------- fill faces (layer patches) --------------------
if isfield(S,'faces') && ~isempty(S.faces) && isstruct(S.faces)

    % ============================================================
    % Case (B): cycles given as edge lists (can overlap in edges)
    % ============================================================
    if isfield(S.faces,'cycles') && ~isempty(S.faces.cycles) && iscell(S.faces.cycles)

        Cyc = S.faces.cycles;
        F   = numel(Cyc);

        for f = 1:F
            cyc0 = Cyc{f};
            if isempty(cyc0)
                continue
            end

            % reorder edges into a chain if possible
            cyc = cellcomplex_order_cycle_edges(cyc0);

            % If reordering fails because of overlaps / shared edges,
            % fall back to a local polygon extraction from cyc0.             
            if isempty(cyc)                                                  
                v = cellcomplex_cycle_vertices_from_edges(cyc0);             
                if isempty(v), continue, end                                 
            else
                v = [cyc(:,1); cyc(end,2)];
            end

            % build polygon coordinates
            xy = Nodes(v,:);
            if v(end) ~= v(1)
                xy = [xy; xy(1,:)];
            end

            % sign from polygon area (y-down convention)
            fs = cellcomplex_cycle_sign_from_area_yDown(xy);

            % patch on top of previous patches (overlaps show naturally)
            if fs > 0
                patch(xy(:,1), xy(:,2), 'r', 'FaceAlpha', 0.04, 'EdgeColor','none');
            elseif fs < 0
                patch(xy(:,1), xy(:,2), 'b', 'FaceAlpha', 0.04, 'EdgeColor','none');
            end
        end

    % ============================================================
    % Case (A): faces given as vertex lists (already polygons)
    % ============================================================
    elseif isfield(S.faces,'vertices') && ~isempty(S.faces.vertices)

        V = S.faces.vertices;
        F = size(V,1);

        if isfield(S.faces,'sign') && ~isempty(S.faces.sign)
            fs = S.faces.sign(:);
        else
            fs = zeros(F,1);
        end

        for f = 1:F
            if fs(f)==0
                continue
            end

            v  = V(f,:);
            xy = Nodes(v,:);
            xy = [xy; xy(1,:)];

            if fs(f) > 0
                patch(xy(:,1), xy(:,2), 'r', 'FaceAlpha', 0.04, 'EdgeColor','none');
            else
                patch(xy(:,1), xy(:,2), 'b', 'FaceAlpha', 0.04, 'EdgeColor','none');
            end
        end

    % ============================================================
    % Fallback: array-of-face structs
    % ============================================================
    elseif numel(S.faces) > 1

        for f = 1:numel(S.faces)

            if isfield(S.faces(f),'sign') && ~isempty(S.faces(f).sign)
                fs = S.faces(f).sign;
            else
                fs = 0;
            end
            if fs==0
                continue
            end

            if isfield(S.faces(f),'verts') && ~isempty(S.faces(f).verts)
                v = S.faces(f).verts(:);
            elseif isfield(S.faces(f),'vertices') && ~isempty(S.faces(f).vertices)
                v = S.faces(f).vertices(:);
            else
                continue
            end

            xy = Nodes(v,:);
            xy = [xy; xy(1,:)];

            if fs > 0
                patch(xy(:,1), xy(:,2), 'r', 'FaceAlpha', 0.04, 'EdgeColor','none');
            else
                patch(xy(:,1), xy(:,2), 'b', 'FaceAlpha', 0.04, 'EdgeColor','none');
            end
        end
    end
end

axis equal;
box off; axis off;                                                  
hold off;

end

% =====================================================================
function cyc = cellcomplex_order_cycle_edges(cyc0)
% Try to reorder edges into a single directed cycle (unique successor).
% Returns [] if it cannot form a single cycle chain (e.g., overlaps/branches).

if isempty(cyc0)
    cyc = cyc0;
    return
end

K = max(cyc0(:));
next = zeros(K,1);

% If multiple outgoing edges from a node exist in cyc0, ordering is ambiguous.
% We detect it and fail fast to allow fallback.                             % <— FIX
for r = 1:size(cyc0,1)                                                   % <— FIX
    u = cyc0(r,1); v = cyc0(r,2);                                         % <— FIX
    if next(u) ~= 0 && next(u) ~= v                                       % <— FIX
        cyc = [];                                                         % <— FIX
        return                                                            % <— FIX
    end                                                                   % <— FIX
    next(u) = v;                                                          % <— FIX
end                                                                       % <— FIX

vstart = cyc0(1,1);
vcur   = vstart;
cyc    = zeros(0,2);

while true
    vnext = next(vcur);
    if vnext==0
        cyc = zeros(0,2);
        break
    end
    cyc(end+1,:) = [vcur vnext]; %#ok<AGROW>
    vcur = vnext;
    if vcur==vstart
        break
    end
end

end

% =====================================================================
function v = cellcomplex_cycle_vertices_from_edges(cyc0)
% Extract a polygon vertex order from an edge set even when the directed
% chain is not globally consistent (e.g., two cycles share an edge).
%
% Strategy:
%   Treat cyc0 as an undirected adjacency, then walk a simple cycle
%   using a greedy "do not immediately go back" rule.
%
% This is intended only for visualization.                                  % <— FIX

if isempty(cyc0)
    v = [];
    return
end

% build undirected adjacency list
K = max(cyc0(:));
Adj = cell(K,1);
for r = 1:size(cyc0,1)
    a = cyc0(r,1); b = cyc0(r,2);
    Adj{a} = [Adj{a} b];
    Adj{b} = [Adj{b} a];
end

% start from first edge
v0 = cyc0(1,1);
v1 = cyc0(1,2);

v = v0;
prev = v0;
cur  = v1;

% walk until we return or fail
maxSteps = 5*size(cyc0,1);
for s = 1:maxSteps
    v(end+1,1) = cur; %#ok<AGROW>
    nbrs = Adj{cur};
    if isempty(nbrs)
        v = [];
        return
    end

    % choose next that is not the immediate predecessor (if possible)
    if numel(nbrs)==1
        nxt = nbrs(1);
    else
        cand = nbrs(nbrs ~= prev);
        if isempty(cand)
            nxt = nbrs(1);
        else
            nxt = cand(1);
        end
    end

    prev = cur;
    cur  = nxt;

    if cur == v0
        break
    end
end

% must close
if v(end) ~= v0
    % try to close by appending v0 if adjacent
    if any(Adj{v(end)} == v0)
        v(end+1,1) = v0;
    else
        v = [];
    end
end

end

% =====================================================================
function sgn = cellcomplex_cycle_sign_from_area_yDown(xy)
x = xy(:,1);
y = xy(:,2);

if x(1)~=x(end) || y(1)~=y(end)
    x = [x; x(1)];
    y = [y; y(1)];
end

A = 0;
for i = 1:(numel(x)-1)
    A = A + x(i)*y(i+1) - x(i+1)*y(i);
end
A = 0.5*A;

% y-axis points down in image coordinates
if A < 0
    sgn = 1;
else
    sgn = -1;
end

end

