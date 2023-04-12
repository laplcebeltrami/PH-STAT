function C = randcorr(p)
% RANDCORR Generate a random p x p correlation matrix
%   C = randcorr(p) Implements algorithm [1] for generating a random
%   p x p correlation matrix using the hyperspherical coordinate representation.
%   The random sampling of angles from sin^k(theta) is based on the algorithm in [2].
%
%   The input argument:
%       p        - size of the correlation matrix
%
%   Return values:
%       C        - [p x p] correlation matrix
%
%   References:
%   [1] Mohsen Pourahmadi and Xiao Wang,
%       Distribution of random correlation matrices: Hyperspherical parameterization of the Cholesky factor,
%       Statistics & Probability Letters, Volume 106, November 2015, Pages 5-12
%
%   [2] Enes Makalic and Daniel F. Schmidt,
%       An efficient algorithm for sampling from $\sin^k(x)$ for generating random correlation matrices
%       arxiv, 2018
% 
%   (c) Copyright Enes Makalic and Daniel F. Schmidt, 2018

%% Error checking
if(numel(p) > 1)
    error('Input p must be a scalar integer greater than 1.');
end
if(isinf(p) || isnan(p))
    error('Input p must be a scalar integer greater than 1.');
end   
if(p < 2)
    error('Input p must be a scalar integer greater than 1.');
end
if(mod(p,1))
    error('Input p must be a scalar integer greater than 1.');
end

%% Step 1 - generate angles theta from PDF (sin(theta))^k, k>=1, 0<theta<pi
e = ones(p,1);
theta = zeros(p,p);
for j = 1:(p-1)
    theta((j+1):p,j) = draw_theta( (p-j)*e((j+1):p) );
end

%% Step 2 - construct lower triangular Cholesky factor
L = ones(p,p);
for i = 2:p
    L(i,2:i) = cumprod( sin(theta(i,1:i-1)) );
end
L = L .* tril(cos(theta));

%% Form correlation matrix
C = L*L';

%% Ensure all diagonal entries are set to 1
C(1:(p+1):end) = 1;

end

%% Draw theta ~ (sin(x))^k, where 0 < x < pi using rejection sampling
function x = draw_theta(k)

k = k(:);
N = length(k);
logconst = 2*log(pi/2);

%% Sampling loop - vectorized
x = zeros(N,1);
accept = false(N,1);
while(~all(accept))
    
    %% index of sampels that need to be accepted 
    ix = ~accept;
    T = sum(ix);

    %% Beta(k+1,k+1) rng
    g1 = randg(k(ix)+1);
    g2 = randg(k(ix)+1);
    x(ix) = pi*(g1 ./ (g1 + g2));         
    
    %% Check acceptance
    accept(ix) = log(rand(T,1))./k(ix) < logconst + log(sin(x(ix))) - log(x(ix)) - log(pi-x(ix));
       
end

end