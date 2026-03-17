function h = plot_graph(G)
% plot_graph(G)
%
% Plot graph in 3D using automatic layout.
%
% INPUT
%   G : MATLAB graph object
%
% OUTPUT
%   h : figure handle
%
% (C) 2026 Moo K. Chung

h = figure;

% --- generate coordinates from graph layout 
p = plot(G);
X = p.XData;
Y = p.YData;
Z = zeros(size(X));   % embed into 3D  
delete(p);            % remove default plot  

% --- plot nodes ---
scatter3(X, Y, Z, 120, 'k', 'filled');

hold on;

% --- node labels ---   
for i = 1:length(X)
    text(X(i)+0.1, Y(i)+0.1, Z(i)+0.1, ...
        ['  ' num2str(i)], ...   
        'FontSize', 18, ...
        'Color', 'k');
end

% --- draw edges ---
E = G.Edges.EndNodes;
for i = 1:size(E,1)
    plot3([X(E(i,1)) X(E(i,2))], ...
          [Y(E(i,1)) Y(E(i,2))], ...
          [Z(E(i,1)) Z(E(i,2))], ...
          'k', 'LineWidth', 3);
end

axis equal;
%set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);   
grid on; axis on

set(gca,'FontSize',12);
set(gcf,'Color','w','InvertHardcopy','off');

end