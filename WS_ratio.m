function observed = WS_ratio(lossMtx, nGroup1, nGroup2)
%function observed = WS_ratio(lossMtx, nGroup1, nGroup2)   
%
% Compute the ratio statistic bwetween-group over within-group given in
%
% [1] Songdechakraiwut, T., Shen, L., Chung, M.K. 2021 Topological learning and 
%its application to multimodal brain network integration, Medical Image 
%Computing and Computer Assisted Intervention (MICCAI), LNCS 12902:166-176 
%
% [2] Songdechakraiwut, T. Chung, M.K. 2023 Topological learning for brain networks, 
% Annals of Applied Statistics 17:403-433, arXiv:2012.00675.
% https://arxiv.org/pdf/2012.00675.pdf
%
%
% INPUT
% lossMtx          : matrix whose entries are pair-wise losses/distances
% nGroup1, nGroup2 : sample size in each group
%
%
% OUTPUT
% observed : observed p-value
%
% %If you are using any part of the code, please reference the above paper.
%
% (C) 2020 Tananun Songdechakraiwut, Moo K. Chung
%          University of Wisconsin-Madison
%
%  Contact tananun@cs.wisc.edu or mkchung@wisc.edu
%  for support/permission with the codes 
%
% Update history
%     2020, November 11 created - Songdechakraiwut
%     2020, December 10 edited - Chung
%     2023, Febuary 10 dedeled last argument that inputs the form of ratio statistic Chung
%
%
% within groups
within = 0;
% sum of pair-wise distances within groups
for i = 1:nGroup1 % group 1
    for j = i + 1:nGroup1
        within = within + lossMtx(i, j);
    end
end
for i = nGroup1 + 1:nGroup1 + nGroup2 % group 2
    for j = i + 1:nGroup1 + nGroup2
        within = within + lossMtx(i, j);
    end
end

% between groups
% sum of pair-wise distances between the groups.

between = 0;
for i = 1:nGroup1
    for j = nGroup1 + 1:nGroup1 + nGroup2
        between = between + lossMtx(i, j);
    end
end

% denomGroup1 = (nGroup1*(nGroup1-1) + nGroup2*(nGroup2-1))/2;
% denomGroup2 = nGroup1*nGroup2;
% Papers [1] and [2] normalize the ratio statistic by the number of
% elements. This is not exaclty needded. We no longer normalize. 

observed = between/within; 

end



