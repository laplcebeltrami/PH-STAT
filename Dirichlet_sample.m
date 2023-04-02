function sample = Dirichlet_sample(alpha,weight,n)
%function sample = Dirichlet_sample(alpha,weight,n)
%
% The function generates samples from a mixture Dirichlet distribution on a 2-simplex
%
% Input:
%    alpha :  3*k matrix, where k is the number of components
%   weight :  k*1 or 1*k vector. To sample from a single dirichlet 
%             distribution, put weight = 1 (a scalar)
%        n :  scalar, sample size
%
% Output:
%   sample :  n*2 matrix. Each row of the output corresponds to one sample.
%
%
% (C) 2022 Zijian Chen, Moo K. Chung
%     University of Wisconsin-Madison
%
%
% To display contour plot for the theoretical density, use function:
% theoreticalDensity.m


weight_sum = cumsum(weight);
sample = zeros(n,2);


for i = 1:n
    U = rand(1);
    index = find( weight_sum > U, 1 );
    sample(i,:) = DirichletComponent(alpha(:,index),1);
end



end


function s3 = DirichletComponent(alpha,n)

% this function generates samples from a single Dirichlet component on a 2-simplex
% (upper triangle in the plane)
% 
% Input: 
%           alpha: 3*1 column vector
%               n: number of samples
% Output:
%              s3: Dirichlet samples
%
% Example:
%       s3 = Dirichlet_2d_sample([3,2,3]',100);




s1 = gamrnd(repmat(alpha',n,1),1,n,3); 
s1 = s1./repmat(sum(s1,2),1,3);

s2 = s1(:,[1,2]);

s3 = [s2(:,1),1-s2(:,2)];

end
