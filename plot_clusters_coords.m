function plot_clusters_coords(coord_d)
% plot_clusters_coords(coord_d)
%
% Plot coordinate observations for each group.
%
% INPUT
% coord_d : cell array (#networks x #clusters) containing coordinates
%
% (C) 2026 Moo K. Chung
% University of Wisconsin-Madison

[nN, nC] = size(coord_d);   % extract dimensions

figure;

nr = ceil(sqrt(nC));        % automatic subplot layout
nc = ceil(nC/nr);

for j = 1:nC
    subplot(nr,nc,j)
    hold on
    
    for i = 1:nN
        obs = coord_d{i,j};
        plot(obs(:,1), obs(:,2), '.k');
    end
    figure_bigger(16);
    axis square
    xlim([-3 3])
    ylim([-3 3])
    title(['Group ' int2str(j)])
end



end