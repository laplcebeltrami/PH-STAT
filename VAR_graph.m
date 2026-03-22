function Z = VAR_graph(Edges, X, p, noise_level)
% VAR_graph generates VAR-based node time series Z whose cross-coupling is
% induced by a directed edge-flow X on the graph Edges.
%
%   Z = VAR_graph(Edges, X, p, noise_level)
%
% Inputs
%   Edges       : m x 2 directed edges [i j] (i -> j), node labels are positive integers
%   X           : m x 1 edge weights aligned to Edges (same ordering)
%   p           : VAR order (positive integer)
%   noise_level : innovation std (scalar)
%
% Output
%   Z           : n x T node time series, where
%                 T = 200 is the simulated length,
%                 n = max(Edges(:)) is the largest node index appearing in Edges.  % <— FIX: explicit definition
%
% Notes
%   - n = max(Edges(:)) assumes nodes are labeled 1..n and that every node index
%   up to n is meaningful. If there are isolated (singleton) nodes that never
%   appear in Edges, they are not represented here: the VAR is simulated only
%   on the edge-induced node set (nodes that appear in Edges).              % <— FIX: clarified scope
%
%   - Lag matrices are built directly from the directed flow matrix W.
%
%   - W(i,j) stores total weight for edge i -> j (accumulated if repeated).
%
%   - Incoming edges drive each node via transpose W': (i->j) contributes to row j.
%
%   - No burn-in, no spectral-radius normalization.
%
% (c) 2025 Moo K. Chung
% University of Wisconsin-Madison
% mkchung@wisc.edu 
% 
% The code is downloaded from 
% https://github.com/laplcebeltrami/hodge
% If you are using the code, refernce one of Hodge papers listed in GitHub.  
%
% Update history: 2025 November 7 created
%                 2026 Feb 7 comments    

% ----- fixed sizes / coefficients -----
T    = 200;        % length of series

% ----- directed flow matrix W: W(i,j) = weight for i -> j -----
n = max(Edges(:)); % number of nodes implied by the largest node index   % <— FIX: comment
W = zeros(n);
for k = 1:size(Edges,1)
    i = Edges(k,1); j = Edges(k,2);
    W(i,j) = W(i,j) + X(k);
end
W(1:n+1:end) = 0;   % no self-edges

% ----- build lag matrices A{1..p} from W (simple, uniform) -----
A = cell(1,p);

rho   = 0.1;
coup  = 0.1;
decay = 0.5;

for ell = 1:p
    A{ell} = (ell==1)*rho*eye(n) + (coup * decay^(ell-1)) * W.';  % incoming drives node
end

% ----- simulate VAR(p): Z_t = sum_{ell=1}^p A{ell} Z_{t-ell} + eps_t -----
Z = zeros(n, T);
Z(:,1:p) = noise_level * randn(n, p);  % initialize with noise (no burn-in)
for t = p+1:T
    acc = zeros(n,1);
    for ell = 1:p
        acc = acc + A{ell} * Z(:, t-ell);
    end
    Z(:, t) = acc + noise_level * randn(n,1);
end
end

% Previous Version
% 
% function Z = VAR_graph(Edges, X, p, noise_level)
% % VAR_graph generates VAR-based directed edge flow along Edges of a graph
% % using edge flow X as the ground truth. 
% %
% %   [Z, A_true] = VAR_graph(Edges, X, p, noise_level)
% %
% % Inputs
% %   Edges       : m x 2 directed edges [i j] (i -> j), nodes labeled 1..n
% %   X           : m x 1 edge weights aligned to Edges
% %   p           : VAR order (positive integer)
% %   noise_level : innovation std (scalar)
% %
% % Outputs
% %%   Z      : n x T time series (fixed at T=200; first p samples are i.i.d. noise)
% %
% %   n = max(Edges(:)) treats node labels as 1..n and sets n to the largest
% %   index appearing in Edges. This assumes every node participates in at least
% %   one edge. If isolated (singleton) nodes exist, they will be excluded and
% %   the simulated VAR will be defined only on the edge-induced node set.  
% %
% %%   A_true : struct with fields:
% %              .A      = 1 x p cell of lag matrices A{1},...,A{p}
% %              .rho    = self-retention coefficient used
% %              .gammas = 1 x p vector of coupling scales per lag
% %
% % Notes
% %   - No burn-in, no spectral-radius normalization.
% %   - Lag matrices are built directly from the directed flow matrix W.
% %   - Incoming edges drive each node via transpose (W'): (i->j) contributes to row j.
% %
% %
% % (c) 2025 Moo K. Chung
% % University of Wisconsin-Madison
% % mkchung@wisc.edu 
% % 
% % The code is downloaded from 
% % https://github.com/laplcebeltrami/hodge
% % If you are using the code, refernce one of Hodge papers listed in GitHub.  
% %
% % Update history: November 7, 2025
% 
% 
%     % ----- fixed sizes / coefficients -----
%     T    = 200;        % length of series
% 
%     % ----- directed flow matrix W: W(i,j) = weight for i -> j -----
%     n = max(Edges(:));
%     W = zeros(n);
%     for k = 1:size(Edges,1)
%         i = Edges(k,1); j = Edges(k,2);
%         W(i,j) = W(i,j) + X(k);
%     end
%     W(1:n+1:end) = 0;   % no self-edges
% 
% 
%     % ----- build lag matrices A{1..p} from W (simple, uniform) -----
%     A = cell(1,p);
% 
%     % example: target ~0.90 with P=4
%     rho   = 0.1;
%     coup  = 0.1;
%     decay = 0.5;
% 
%     %rho   small self-retention (keeps series from being too smooth)
%     %coup   overall cross-coupling scale (smaller -> harder task)
%     % Larger coup: Stronger influence of flows between nodes. Each node’s next state 
%     % depends more heavily on its neighbors’ past states. This amplifies the causal 
%     % structure, makes cycles more coherent, and usually boosts performance 
%     % (cosine similarity closer to 1). But if coup is too large, the dynamics 
%     % can become overly smooth or even unstable.
% 	% Smaller coup: Weaker cross-coupling. The process is more noise-driven, 
%     % causal patterns are harder to detect, and performance drops.
% 
%     %decay  geometric decay per lag (0<decay<1)
%     %large decay : long memory, smoother series, higher cosine accuracy,
%     %small decay : short memory, rougher series, lower accuracy.
% 
%     for ell = 1:p
%         A{ell} = (ell==1)*rho*eye(n) + (coup * decay^(ell-1)) * W.';  % fill every A{ell}
%     end
% 
%     % ----- simulate VAR(p): Z_t = sum_{ell=1}^p A{ell} Z_{t-ell} + eps_t -----
%     Z = zeros(n, T);
%     Z(:,1:p) = noise_level * randn(n, p);  % initialize with noise (no burn-in)
%     for t = p+1:T
%         acc = zeros(n,1);
%         for ell = 1:p
%             acc = acc + A{ell} * Z(:, t-ell);
%         end
%         Z(:, t) = acc + noise_level * randn(n,1);
%     end
% 
% 
% end
% 
