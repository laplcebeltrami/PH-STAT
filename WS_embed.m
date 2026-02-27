function Z = WS_embed(C)
% function Z = WS_embed(C)
%
% Performs topological embedding through the birth-death decomposition.
%
% INPUT
%   C : connectivity matrices (#nodes x #nodes x #subjects)
%
% OUTPUT
%   Z : structured variable containing
%       Z.x      - x-coordinates (birth embedding)
%       Z.y      - y-coordinates (death embedding)
%       Z.center - embedding center
%       Z.mean_b - mean birth values
%       Z.mean_d - mean death values
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
%
% (C) 2023 Moo K. Chung
%     University of Wisconsin-Madison
%     mkchung@wisc.edu
%
% 2026 Feb 27, ouput structured

% Birth-death decomposition
[Wb, Wd] = WS_decompose(C);

birth = squeeze(Wb(:,3,:));
death = squeeze(Wd(:,3,:));

n_sub = size(birth,2);


% Centering
mean_b = mean(birth,2);
mean_d = mean(death,2);

mean_bmatrix = repmat(mean_b, 1, n_sub);
mean_dmatrix = repmat(mean_d, 1, n_sub);


% Embedding coordinates
Z.x = mean(birth - mean_bmatrix, 1);
Z.y = mean(death - mean_dmatrix, 1);


% Embedding center
Z.center = [mean(mean_b), mean(mean_d)];

% Optional: store additional useful quantities
Z.mean_b = mean_b;
Z.mean_d = mean_d;

end

% function [z, c] = WS_embed(C);
% %function [z c] = WS_embed(C);
% %    
% % Performs topological embeeding through the birth-death decomposition
% %
% % INPUT
% %   C         : collection of connectivity matrice of size 
% %               # of nodes x # of nodes x # subjects
% %
% % OUTPUT
% % z.x and z.y : x- and y-coordiantes of the topological embdding. 
% % c           : center of cembedding. 
% %
% % The method is published in
% % 
% % [1] Songdechakraiwut, T., Shen, L., Chung, M.K. 2021 Topological learning and 
% %its application to multimodal brain network integration, Medical Image 
% %Computing and Computer Assisted Intervention (MICCAI), LNCS 12902:166-176 
% %
% % [2] Songdechakraiwut, T. Chung, M.K. 2020 Topological learning for brain 
% % networks, arXiv: 2012.00675. 
% % 
% %
% % If you are using any part of the code, please reference the above paper.
% % The function is downloaded from 
% % http://pages.stat.wisc.edu/~mchung/publication.html
% %
% %
% % (C) 2023 Moo K. Chung
% %     University of Wisconsin-Madison
% %  Contact mkchung@wisc.edu for support 
% %
% % Update history
% %   2023 April 16
% 
% %Birth-death decomposition
% [Wb Wd] = WS_decompose(C);
% % Wb      : birth edge set  (p-1) x 3 x n, where p is # of nodes and n is # of subjects
% % Wd      : death edge set  (p-1)*(p-2)/2 x 3 x n, where p is # of nodes and n is # of subjects
% 
% 
% birth = squeeze(Wb(:,3,:)); %birth values
% death = squeeze(Wd(:,3,:)); %death values
% 
% n_sub= size(birth,2); %number of subjects
% 
% mean_b = mean(birth,2); %mean of birth values over all subjects
% mean_d = mean(death,2); %mean of death values over all subjects
% 
% mean_bmatrix = repmat(mean_b, 1, n_sub);
% mean_dmatrix = repmat(mean_d, 1, n_sub); 
% 
% z.x = mean(birth-mean_bmatrix,1);
% z.y = mean(death-mean_dmatrix,1);
% 
% c=[mean(mean_b) mean(mean_d)];
% 
% 
