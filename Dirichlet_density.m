function Z = Dirichlet_density(alpha, weight, x)
%function Z = Dirichlet_density(x)
%
% Given the coordinates in the persistent diagram, the function computes the theoretical 
% Dirichlet distribution with parameters alpha and mixing weight on the persistent diagram. 
%
% INPUT:  
%    alpha :  3*k matrix, where k is the number of components
%   weight :  k*1 or 1*k vector. To sample from a single dirichlet 
%             distribution, put weight = 1 (a scalar)
%       x   coordinates. x(:,1) is the x-coordinates, y(:,2)is the y-coordinates
%
%
% The code is downloaded from
% https://github.com/laplcebeltrami/PH-STAT
%
%
% (C) 2022 Zijian Chen, Moo K. Chung
%     University of Wisconsin-Madison
%
%
% Modified from Dirichelt_theoretical.m

alpha=alpha';
D = inline('gamma(a1+a2+a3)./(gamma(a1)*gamma(a2)*gamma(a3)).*x1.^(a1-1).*(1-x2).^(a2-1).*(x2-x1).^(a3-1)');


X=x(:,1);
Y=x(:,2);

% compute density
Z = zeros(size(X));
for i = 1:length(weight)
    Z = Z + weight(i)*D(alpha(i,1),alpha(i,2),alpha(i,3),X,Y);
end

% for i = 1:Nplot
%     for j = Nplot:-1:i
%         Z(i,j)=0;
%     end
% end
% draw contour
%contourf(X,Y,Z,100,'LineColor','none')


end

