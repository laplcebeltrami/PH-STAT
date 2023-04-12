function pvalue = online_pvalue(stat, observed)
%function pvalue = online_pvalue(stat, observed)
%
% Computes the pvalue based on the collection of statistics value 
% and observation sequentically. It is needed for various applications
% where we need to know how p-values change. 
% 
% INPUT
%       stat:     sequence of test statistic of size n (# of statistics) x l (# of variables/connections)
%       observed: observed test statistic
%
% OUTPUT
%       pvalue  : sequence of p-values    
% 
% 
% The seqencial p-value 
% computation algorithm is given in 
%
% [1] Chung, M.K., Xie, L., Huang, S.-G., Wang, Y., Yan, J., Shen, L. 2019. 
% Rapid Acceleration of the permutation test via transpositions, International 
% Workshop on Connectomics in NeuroImaging, Lecture Notes in Computer Science 11848:42-53. 
% A different version arXiv:1812.06696
%
% [2] Songdechakraiwut, T. Chung, M.K. 2020 Topological learning for brain networks, arXiv: 2012.00675.
%
% This code is downloaded from
% http://www.stat.wisc.edu/~mchung/transpositions
%
% (C) 2019 Moo K. Chung  mkchung@wisc.edu
% University of Wisconsin-Madison
%
% Update histroy: 2019 May created
%                 2021 Oct. 11 Validation done 
% 
%
% Validaiton against built-in Matlab function normrnd.m with 10000 N(0,1)
% random numbers. 
% norminv(0.95,0,1) = 1.6449 treshold corresponding to pvalue = 1-0.95 probability
% norminv(0.5,0,1) = 0
% norminv(0.05,0,1) = -1.6449
%
% Based on the above ground truth, we should get p-values
%
% stat=normrnd(0,1,10000,1);
% observed =-1.6449  %1.6449  
% pvalue = online_pvalue(stat, observed)
% pvalue(end) 
% 0.0482  0.4889   0.0489


% pvalue can be computed iteratively as
% pvalue(i+1) = (pvalue(i) * i  + ( t(i+1)>=observed ) / (i+1)


n=size(stat,1); % number of statistics.
l=size(stat,2); % numver of variables/connections

pvalue=zeros(n,l);

if observed>=mean(stat) %right tail
    
    pvalue(1,:)= (stat(1)>=observed); %initial p-value. It is either 0 or 1.
    for i=2:n
        pvalue(i,:) = (pvalue(i-1,:) * (i-1) + (stat(i,:)>=observed))/i;
    end
    
else %left tail
    
    pvalue(1,:)= (stat(1)<=observed); %initial p-value. It is either 0 or 1.
    for i=2:n
        pvalue(i,:) = (pvalue(i-1,:) * (i-1) + (stat(i,:)<=observed))/i;
    end
    
end

%convergence plot for p-value
figure; 
set(gcf, 'Position', [400 400 600 250])


plot(pvalue(1:10000),'k','LineWidth', 2)
xlabel ('Number of permutations'); ylabel('p-value')
whitebg(gcf,'w');
set(gcf,'Color','w','InvertHardcopy','off');

set(gca, 'fontsize',16) 
