function Z = Dirichlet_theoretical(alpha,weight)
%function Z = Dirichlet_theoretical(alpha,weight)
%
% The function computes the theoretical Dirichlet distribution on the persistent diagram. 
%
% (C) 2022 Zijian Chen, Moo K. Chung
%     University of Wisconsin-Madison

alpha=alpha';
D = inline('gamma(a1+a2+a3)./(gamma(a1)*gamma(a2)*gamma(a3)).*x1.^(a1-1).*(1-x2).^(a2-1).*(x2-x1).^(a3-1)');

% create grids
Nplot = 100;
x = linspace(0,1,Nplot);
y = linspace(0,1,Nplot);
[X,Y] = meshgrid(x,y);

% compute density
Z = zeros(size(X));
for i = 1:length(weight)
    Z = Z + weight(i)*D(alpha(i,1),alpha(i,2),alpha(i,3),X,Y);
end

for i = 1:Nplot
    for j = Nplot:-1:i
        Z(i,j)=0;
    end
end

% draw contour
contourf(X,Y,Z,100,'LineColor','none')


end

