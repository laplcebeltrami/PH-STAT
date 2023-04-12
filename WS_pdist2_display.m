function  WS_pdist2_display(con_i, con_j)
%function  WS_pdist2_display(con_i, con_j)
%
% The function displays the Wasserstein distance matrix.          

%
% (C) 2023- Moo K. Chung                         
%      University of Wisconsin-Madison
%      mkchung@wisc.edu
%

si = size(con_i);
con_iresahpe = reshape(con_i, si(1)*si(2),si(3));

sj = size(con_j);
con_jresahpe = reshape(con_j, sj(1)*sj(2),sj(3));

con = [con_iresahpe con_jresahpe]';
L2loss = pdist2(con,con);


lossMtx = WS_pdist2(con_i,con_j);


figure; 
subplot(2,2,1); imagesc(L2loss); colorbar; title('Euclidean dist.')
set(gca,'FontSize',12); 
axis square


subplot(2,2,2); imagesc(lossMtx.D0); colorbar; title('Wasserstein dist. 0D')
set(gca,'FontSize',12); 
axis square

subplot(2,2,3); imagesc(lossMtx.D1); colorbar; title('Wasserstein dist. 1D')
set(gca,'FontSize',12); 
axis square

subplot(2,2,4); imagesc(lossMtx.D01); colorbar;  title('Wasserstein dist. 0D+1D')
set(gca,'FontSize',12); 
axis square



set(gcf,'Color','w','InvertHardcopy','off');