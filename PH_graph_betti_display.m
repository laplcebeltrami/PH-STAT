function  PH_graph_betti_display(con_i, con_j,thresholds)
%function PH_graph_betti_display(beta, threholds)
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

nGroup_i = size(con_i,3);
nGroup_j = size(con_j,3);

betti0_i=[]; betti1_i=[];
for i=1:nGroup_i
    beta = PH_betti(con_i(:,:,i), thresholds);
    betti0_i=[betti0_i; beta.zero];
    betti1_i=[betti1_i; beta.one];
end
mean0_con_i = mean(betti0_i,1); %average Betti-0 curve for i-th group
mean1_con_i = mean(betti1_i,1); %average Betti-1 curve for i-th group

betti0_j=[]; betti1_j=[];
for j=1:nGroup_j
    beta = PH_betti(con_j(:,:,j), thresholds);
    betti0_j=[betti0_j; beta.zero];
    betti1_j=[betti1_j; beta.one];
end
mean0_con_j = mean(betti0_j,1); %average Betti-0 curve for j-th group
mean1_con_j = mean(betti1_j,1); %average Betti-1 curve for j-th group

% visulization of Betti curves. The curves are extremly stable under
% perturbation demonstarting the robustness of the graph filration. Such
% robustnes will not be obtained if you use different filtrations such as
% Rips. 

figure; 
subplot(2,1,1); %Betti-0
for i=1:nGroup_i
     hold on; plot(thresholds,betti0_i(i,:), 'Color', [1 0.6 0.2] , 'LineWidth',1,'LineStyle','-.')
end
hold on; plot(thresholds,mean0_con_i, 'r', 'LineWidth',2,'LineStyle','-')

for j=1:nGroup_j
     hold on; plot(thresholds,betti0_j(j,:), 'Color', [0.3 0.7 1] , 'LineWidth',1,'LineStyle','-.')
end
hold on; plot(thresholds,mean0_con_j, 'b', 'LineWidth',2,'LineStyle','-')

set(gca,'FontSize',18);
xlabel('Correlations');
ylabel('\beta_0'); box on;


subplot(2,1,2); %Betti-1
for i=1:nGroup_i
     hold on; plot(thresholds,betti1_i(i,:), 'Color', [1 0.6 0.2] , 'LineWidth',1,'LineStyle','-.')
end
hold on; plot(thresholds,mean1_con_i, 'r', 'LineWidth',2,'LineStyle','-')

for j=1:nGroup_j
     hold on; plot(thresholds,betti1_j(j,:), 'Color', [0.3 0.7 1] , 'LineWidth',1,'LineStyle','-.')
end
hold on; plot(thresholds,mean1_con_j, 'b', 'LineWidth',2,'LineStyle','-')

set(gca,'FontSize',18);
xlabel('Correlations');
ylabel('\beta_1'); box on;

set(gcf,'Color','w','InvertHardcopy','off');