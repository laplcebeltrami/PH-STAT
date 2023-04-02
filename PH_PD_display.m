function D = PH_PD_display(pairs)
%function D = PH_PD_display(pairs)
%
% Displays persistent diagrms from (birth, death) pairs
%
%
% (C) 2023 Moo K. Chung, University of Wisconsin-Madison
%     Email: mkchung@wisc.edu
%
% The code is downloaded from
% https://github.com/laplcebeltrami/PH-STAT


figure; plot(pairs(:,1),pairs(:,2), 'o','MarkerEdgeColor','k', 'MarkerFaceColor',...
            [0.7 0.7 0.7],'MarkerSize',7)
xlabel('Births'); ylabel('Deaths')

set(gca, 'Fontsize',16);
whitebg(gcf,'w');
set(gcf,'Color','w','InvertHardcopy','off');
