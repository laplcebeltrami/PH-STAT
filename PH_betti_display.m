function  PH_betti_display(beta,thresholds)
%function PH_betti_display(beta)
%
% The function computes displays beta-curves given by beta.zero (betti-0)
% and beta.one (betti-1) over thresholds
%
%          
% The matheamtical details of the methods are published in [1] and [2]. 
% If you are using this code, please reference [1] or [2]. 
%
%% [1] Chung, M.K., Lee, H. Ombao. H., Solo, V. 2019 Exact topological inference 
%%     of the resting-state brain networks in twins, Network Neuroscience 3:674-694
%%     http://www.stat.wisc.edu/~mchung/papers/chung.2019.NN.pdf
%
% [2] Chung, M.K., Huang, S.-G., Gritsenko, A., Shen, L., Lee, H. 2o19
%    Statistical inference on the number of cycles in brain networks. 
%    IEEE International Symposium on Biomedical Imaging (ISBI) 113-116 
%
% Given two Betti curves, we can perform the Exact Topologica Inference (ETI) that provide
% statistical signficance (p-value) of testing the equivalencce of two
% curves. ETI is introduced in
%
% [3] Chung, M.K., Vilalta, V.G., Lee, H., Rathouz, P.J., Lahey, B.B., Zald, D.H. 
%     2017 Exact topological inference for paired brain networks via persistent 
%     homology. Information Processing in Medical Imaging (IPMI) 10265:299-310
%     http://www.stat.wisc.edu/~mchung/papers/chung.2017.IPMI.pdf
%
% [4] Chung, M.K., Luo, Z., Leow, A.D., Alexander, A.L., Richard, D.J., Goldsmith, H.H. 
%     2018 Exact Combinatorial Inference for Brain Images, Medical Image Computing and 
%     Computer Assisted Intervention (MICCAI), 11070:629-637
%     http://www.stat.wisc.edu/~mchung/papers/chung.2018.MICCAI.pdf
%
% [5] Chung, M.K. Lee, H., Gritsenko, A., DiChristofano, A., Pluta, D. 
%     Ombao, H. Solo, V. Topological Brain Network Distances, ArXiv 1809.03878
%     http://arxiv.org/abs/1809.03878
%
%
% (C) 2023- Moo K. Chung                         
%      University of Wisconsin-Madison
%      mkchung@wisc.edu
%

figure;
subplot(2,1,1); plot(thresholds,beta.zero, 'k' , 'LineWidth',2,'LineStyle','-')
title('Betti-0 curve')
set(gca, 'Fontsize',18);
set(gcf,'Color','w','InvertHardcopy','off');
hold on; plot(thresholds,beta.zero, '.k', 'MarkerFaceColor','k',...
                       'MarkerSize',15)


subplot(2,1,2); plot(thresholds,beta.one, 'k', 'LineWidth',2,'LineStyle','-')
title('Betti-1 curve')
set(gca, 'Fontsize',18);
set(gcf,'Color','w','InvertHardcopy','off');
hold on; plot(thresholds,beta.one, '.k','MarkerFaceColor','k',...
                       'MarkerSize',15)