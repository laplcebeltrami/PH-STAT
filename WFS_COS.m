function [wfs beta]=WFS(Y,x,k,sigma)
%
% FUNCTION [wfs beta]=WFS(Y,x,k,sigma)
%
% Construct the k-th degree weighted sine and cosine series 
% representation on functions Y(x). This is equivalent 
% to performing heat kernel smoothing with bandwidth
% sigma. 
%
% Y    : vector of functions (n x m) n: number of sampling points, m: number of input functions. 
% x    : domain of function Y. it has to be normalized to [0,1]. size(Y,1) must qual to the lenght of X.
% k    : degree of expansion
% sigma: amount of smoothing
%
% The technical detail and explanation on paramters are given in [1] and [2]
% 
% [1] Chung, M.K., Adluru, N., Lee, J.E., Lazar, M.,Lainhart, J.E., Alexander, A.L.
%     2010. Cosine series representation of 3D curves and its application to white 
%     matter fiber bundles in diffusion tensor imaging. Statistics and Its Interface. 
%     3:69-80. http://www.stat.wisc.edu/~mchung/papers/chung.2010.SII.pdf
%
% [2] Chung, M.K., Singh, V., Kim, P.T., Dalton, K.M., Davidson, R.J. 2009.
%     Topological characterization of signal in brain images using the min-max diagram.
%     12th International Conference on Medical Image Computing and Computer Assisted 
%     Intervention (MICCAI). Lecture Notes in Computer Science (LNCS). 5762:158-166.
%     http://www.stat.wisc.edu/~mchung/papers/miccai.2009.pdf 
%
% (C) Moo K. Chung 2008
%     mkchung@wisc.edu
%     University of Wisconsin-Madison
%
%  2008 Created
%  2017 Oct. 27 vectorized version, documantation updated
%  2019 June 2. Function name changed from WFS_COS_V2 to WFS_COS
%  2023 Apri 1  Comments updated


%Size of data
m=size(Y,2);

% Cosine basis
psi1=inline('sqrt(2)*cos(pi*l.*x)');

% Sine basis
psi2=inline('sqrt(2)*sin(pi*l.*x)');


% design matrix
DX1=[];
for l=0:k
    psi=psi1(l,x);
    DX1=[DX1 psi];
end

%Include this part for sine basis as well. 
% DX2=[];
% for l=1:k
%    psi=psi2(l,x);
%    DX2=[DX2 psi];
% end
% DX=[DX1 DX2];

DX=DX1;

beta=pinv(DX'*DX)*DX'*Y;

weight = exp(-[0:k]'.^2*pi^2*sigma);
weight = repmat(weight, [1,m]);
%Both sine and cosine basis
%weight = [exp(-[0:k]'.^2*pi^2*sigma) ;exp(-[1:k]'.^2*pi^2*sigma)];
wfs=DX*(weight.*beta);
%wfs = DX*beta;

