function connectivity_stack_display(con)
% connectivity_stack_display  Display connectivity matrices as a stack.
%
% INPUT
%   con : p x p x n array of connectivity matrices
%
% Each connectivity matrix has identical display size, a black boundary,
% and a subject label.
%
% (C) 2026 Moo K. Chung
% University of Wisconsin-Madison

[p,~,n] = size(con);

figure('Color','w');

% Fixed geometry for every matrix
side    = 0.45;     % approximately 2/3 of the previous size 0.68
offset  = 0.025;
left0   = 0.08;
bottom0 = 0.16;

for k = 1:n

    % Shift position only; width and height remain identical
    left   = left0 + offset*(k-1);
    bottom = bottom0 + offset*(n-k);

    ax = axes('Units','normalized', ...
              'Position',[left bottom side side]);

    imagesc(ax,con(:,:,k));
    clim(ax,[-1 1]);

    set(ax, ...
        'XTick',[], ...
        'YTick',[], ...
        'DataAspectRatio',[1 1 1], ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'PositionConstraint','innerposition');

    xlim(ax,[0.5 p+0.5]);
    ylim(ax,[0.5 p+0.5]);

    hold(ax,'on');

    % Black boundary
    plot(ax,[0.5 p+0.5 p+0.5 0.5 0.5], ...
            [0.5 0.5 p+0.5 p+0.5 0.5], ...
            'k','LineWidth',2);

    % Subject label
    text(ax,1.5,1.5,sprintf('Subject %d',k), ...
        'FontSize',11, ...
        'FontWeight','bold', ...
        'Color','k', ...
        'BackgroundColor','w', ...
        'Margin',1, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top');

    hold(ax,'off');
end

colormap(parula);

cb_left   = left0 + side + offset*(n-1) - 0.045;   % move further left
cb_width  = 0.020;
cb_height = 0.67*side;
cb_bottom = bottom0 + (side-cb_height); %/2;

axColor = axes('Units','normalized', ...
               'Position',[cb_left cb_bottom cb_width cb_height], ...
               'Visible','off');

clim(axColor,[-1 1]);
colormap(axColor,parula);

cb = colorbar(axColor);
cb.Position = [cb_left cb_bottom cb_width cb_height];

% Put tick labels on the right
cb.Location = 'eastoutside';
cb.TickDirection = 'out';
cb.AxisLocation = 'out';

cb.FontSize = 14;
cb.Ticks = [-1 -0.5 0 0.5 1];
end
