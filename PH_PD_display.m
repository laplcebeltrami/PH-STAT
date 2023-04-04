function PH_PD_display(pairs, f)
%function D = PH_PD_display(pairs)
%
% Displays persistent diagrms from (birth, death) pairs.
%
% INPUT
%    pairs  : pairs(:,1) is the birth values and pairs(:,2) is the death values
%    f      : functional data on the birth and death values
%
% (C) 2023 Moo K. Chung, University of Wisconsin-Madison
%     Email: mkchung@wisc.edu
%
% The code is downloaded from
% https://github.com/laplcebeltrami/PH-STAT

if nargin <=1

    plot(pairs(:,1),pairs(:,2), 'o','MarkerEdgeColor','k', 'MarkerFaceColor',...
        [0.7 0.7 0.7],'MarkerSize',7)
    xlabel('Births'); ylabel('Deaths')

    set(gca, 'Fontsize',16);
    whitebg(gcf,'w');
    set(gcf,'Color','w','InvertHardcopy','off');

else %if function value f is available
    % change 20 to somehting else to change the size of nodes
    scatter(pairs(:,1), pairs(:,2),20, f,'filled'); colorbar
    axis square; 
    xlabel('Births'); ylabel('Deaths')

    set(gca, 'Fontsize',16);
    whitebg(gcf,'w');
    set(gcf,'Color','w','InvertHardcopy','off');
end
