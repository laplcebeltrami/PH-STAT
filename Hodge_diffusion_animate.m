function Hodge_diffusion_animate(S, X0, Z, dt, nSteps, stride, freqScale)
% Hodge_diffusion_animate
% Animate 1-Hodge diffusion on edge flows:
%   dX(t)/dt + freqScale * Delta1 X(t) = 0
% using closed-form spectral diffusion in the eigenbasis of Delta1.
%
% Inputs:
%   S         : struct with .nodes, .edges, (optional .faces)
%   X0        : initial edge flow [nE x 1]
%   Z         : optional overlay data or [] (kept static)
%   dt        : time step
%   nSteps    : number of simulation steps
%   stride    : plot every stride steps
%   freqScale : multiplies Delta1 to speed diffusion (e.g., 1~20)
%
% This is a minimum revision of animate_edgeflow_EL2_dissipative:
%
% (C) 2026 Moo K. Chung
% University of Wisconsin-Madison
% Update: 2026 Mar 20 created

if nargin < 6
    stride = 1;
end
if nargin < 7
    freqScale = 1;
end

Nodes = S.nodes;
Edges = S.edges;
nE = size(Edges,1);

% build 1-Hodge Laplacian 
B1 = build_B1(size(Nodes,1), Edges);
B2 = build_B2(S, Edges);
if isempty(B2)
    Delta1 = B1' * B1;
else
    Delta1 = B1' * B1 + B2 * B2';
end
Delta1 = freqScale * Delta1;

% eigen-decomposition
[Psi, Lam] = eig(full(Delta1));
lambda = diag(Lam);

% initial condition in modal coordinates 
X = X0(:);
c = Psi' * X;   % Fourier coefficients of initial flow

% color scale 
vmin0 = 0;
vmax0 = max(abs(X0)) * 1.2;
if vmax0 == 0
    vmax0 = 1;
end

% -------------------- draw base once --------------------
figure; clf;
cellcomplex_display(S); hold on;
if ~isempty(Z)
    graph_overlay_timeseries(Nodes, Z);
end
cellcomplex_display_flows_colored(S, X, vmin0, vmax0);
title(sprintf('1-Hodge Diffusion  (t = %.2f)', 0), 'FontSize', 18);
drawnow;
pause(0.03);   

% -------------------- animate --------------------
for step = 1:nSteps

    t = step * dt;

    % exact spectral diffusion update
    X = Psi * (exp(-t * lambda) .* c);

    if mod(step, stride) == 0
        cla;
        cellcomplex_display(S); hold on;

        if ~isempty(Z)
            graph_overlay_timeseries(Nodes, Z);
        end

        cellcomplex_display_flows_colored(S, X, vmin0, vmax0);

        % diagnostic: non-harmonic components decay under diffusion
        r = norm(Delta1 * X); %#ok<NASGU>

        title(sprintf('1-Hodge Diffusion  (t = %.2f)', t), 'FontSize', 18);
        drawnow;
        pause(0.03);   
    end
end

end

% helper function 
% ===============
function B1 = build_B1(nV, Edges)
nE = size(Edges,1);
B1 = sparse(nV, nE);
for e = 1:nE
    tail = Edges(e,1);
    head = Edges(e,2);
    B1(tail,e) = -1;
    B1(head,e) = +1;
end
end

function B2 = build_B2(S, Edges)
B2 = [];
if ~isfield(S,'faces') || isempty(S.faces)
    return
end
F = S.faces;
if isfield(F,'edgeCycles') && ~isempty(F.edgeCycles)
    nF = numel(F.edgeCycles);
    nE = size(Edges,1);
    B2 = sparse(nE, nF);
    for f = 1:nF
        ec = F.edgeCycles{f};
        for r = 1:size(ec,1)
            eIdx = ec(r,1);
            sgn  = ec(r,2);
            B2(eIdx,f) = sgn;
        end
    end
elseif isfield(F,'edges') && ~isempty(F.edges)
    nF = numel(F.edges);
    nE = size(Edges,1);
    B2 = sparse(nE, nF);
    for f = 1:nF
        eList = F.edges{f}(:);
        B2(eList,f) = 1;
    end
end
end