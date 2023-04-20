function [z, c] = WS_embed(C);
%function [z c] = WS_embed(C);
%    
% Performs topological embeeding through the birth-death decomposition
%
% INPUT
%   C         : collection of connectivity matrice of size 
%               # of nodes x # of nodes x # subjects
%
% OUTPUT
% z.x and z.y : x- and y-coordiantes of the topological embdding. 
% c           : center of cembedding. 
%
% The method is published in
% 
% [1] Songdechakraiwut, T., Shen, L., Chung, M.K. 2021 Topological learning and 
%its application to multimodal brain network integration, Medical Image 
%Computing and Computer Assisted Intervention (MICCAI), LNCS 12902:166-176 
%
% [2] Songdechakraiwut, T. Chung, M.K. 2020 Topological learning for brain 
% networks, arXiv: 2012.00675. 
% 
%
% If you are using any part of the code, please reference the above paper.
% The function is downloaded from 
% http://pages.stat.wisc.edu/~mchung/publication.html
%
%
% (C) 2023 Moo K. Chung
%     University of Wisconsin-Madison
%  Contact mkchung@wisc.edu for support 
%
% Update history
%   2023 April 16

%Birth-death decomposition
[Wb Wd] = WS_decompose(C);
% Wb      : birth edge set  (p-1) x 3 x n, where p is # of nodes and n is # of subjects
% Wd      : death edge set  (p-1)*(p-2)/2 x 3 x n, where p is # of nodes and n is # of subjects


birth = squeeze(Wb(:,3,:)); %birth values
death = squeeze(Wd(:,3,:)); %death values

n_sub= size(birth,2); %number of subjects

mean_b = mean(birth,2); %mean of birth values over all subjects
mean_d = mean(death,2); %mean of death values over all subjects

mean_bmatrix = repmat(mean_b, 1, n_sub);
mean_dmatrix = repmat(mean_d, 1, n_sub); 

z.x = mean(birth-mean_bmatrix,1);
z.y = mean(death-mean_dmatrix,1);

c=[mean(mean_b) mean(mean_d)];


