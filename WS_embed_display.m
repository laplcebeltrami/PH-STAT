function WS_embed_display(Z_i, Z_j)
% WS_EMBED_DISPLAY  Visualize two topological embeddings with legend.
%
% INPUT:
%   Z_i, Z_j : structured outputs from WS_embed
%              containing fields:
%              .x, .y, .center
%
% The topological embedding is explained in 
% Chung, M.K., Ramos, C.G., De Paiva, F.B., Mathis, J., Prabhakaran, V., 
% Nair, V.A., Meyerand, M.E., Hermann, B.P., Binder, J.R. and Struck, A.F., 
% 2023. Unified topological inference for brain networks in temporal lobe 
% epilepsy using the Wasserstein distance. NeuroImage, 284, p.120436. arXiv:2302.06673. 
% The simultation setting below comes from 
% Chung, M.K., Huang, S.G., Carroll, I.C., Calhoun, V.D. and Goldsmith, H.H., 
% 2024. Topological state-space estimation of functional human brain networks. 
% PLOS Computational Biology, 20(5), p.e1011869. arXiv:2201:00087
%
% (C) 2025 Moo K. Chung
% University of Wisconsin-Madison
% mkchung@wisc.edu

figure; hold on;

% --- Group i ---
h1 = plot(Z_i.x + Z_i.center(1), ...
          Z_i.y + Z_i.center(2), ...
          'ob', 'MarkerSize',10, 'MarkerFaceColor','g');

h2 = plot(Z_i.center(1), ...
          Z_i.center(2), ...
          'sb', 'MarkerSize',12, 'MarkerFaceColor','b');

% --- Group j ---
h3 = plot(Z_j.x + Z_j.center(1), ...
          Z_j.y + Z_j.center(2), ...
          'or', 'MarkerSize',10, 'MarkerFaceColor','y');

h4 = plot(Z_j.center(1), ...
          Z_j.center(2), ...
          'sr', 'MarkerSize',12, 'MarkerFaceColor','r');

% --- Formatting ---
axis square;
xlabel('Births');
ylabel('Deaths');
title('Topological Embedding');

legend([h1 h2 h3 h4], ...
       {'Group 1', ...
        'Group 1 centroid', ...
        'Group 2', ...
        'Group 2 centroid'}, ...
       'Location','best');

set(gca,'FontSize',16);
set(gcf,'Color','w','InvertHardcopy','off');

if exist('figure_bigger','file')
    figure_bigger(20);
end

end