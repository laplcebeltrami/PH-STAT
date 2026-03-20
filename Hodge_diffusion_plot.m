function hodge_diffusion_plot(A, Phi, lambda, X0)
% hodge_diffusion_plot(A, Phi, lambda, X0)
%
% Plot Hodge diffusion of edge flows using initial flow X0
%
% INPUT
%   A      : adjacency matrix
%   Phi    : eigenvectors of Hodge Laplacian
%   lambda : eigenvalues (vector)
%   X0     : initial edge flow
%
% (C) 2026 Moo K. Chung
% Updates: 2026 Mar 20 created

nE = length(lambda);

% Fourier coefficients
c = Phi' * X0;   % Fourier coefficients

% Time grid
tGrid = linspace(0,4,100);

% Diffusion
Xdiff = zeros(nE, length(tGrid));

for it = 1:length(tGrid)
    t = tGrid(it);
    Xdiff(:,it) = Phi * (exp(-t*lambda) .* c);
end

% Build edge list
[i,j] = find(triu(A,1));
edges = [i j];

% Plot
figure; hold on;

for e = 1:nE
    plot(tGrid, Xdiff(e,:), 'LineWidth', 2);
end

xlabel('Time');
ylabel('Edge flow');
title('Hodge diffusion on edges');

% Edge labels (vertex-based)
labels = arrayfun(@(k) sprintf('(%d,%d)', edges(k,1), edges(k,2)), ...
                  1:nE, 'UniformOutput', false);

legend(labels, 'Location','eastoutside');

set(gca,'FontSize',16);
set(gcf,'Color','w','InvertHardcopy','off');

end