function EC = EC_t_square(u, a, sigma, nu)
%function EC = EC_t_square(u, a, sigma, nu)

% This function computes the expected Euler characteristic for a square domain
% with side length a, given a vector of thresholds u, smoothing sigma, and
% degrees of freedom nu for the T-field.
%
% Inputs:
%   u - A vector of thresholds
%   a - The side length of the square domain
%   sigma - The smoothing parameter
%   nu - The degrees of freedom for the T-field
%
% Output:
%   EC - The expected Euler characteristic for each threshold in u
%
%
% The method and prameters are explained in
% https://arxiv.org/pdf/2405.07835
%
%
% (C) 2024 Moo K. Chung mkchung@wisc.edu
%     Univeristy of Wisconsin-Madison


% Preallocate the EC vector
EC = zeros(size(u));

% Minkowski functionals for square S
mu_0 = 1;
mu_1 = 2 * a;
mu_2 = a^2;

% Loop over each threshold value
for i = 1:length(u)
    % EC densities for the T-field
    rho_0 = 1 - tcdf(u(i), nu);
    rho_1 = (1 / sqrt(2 * sigma^2)) * (1 / (2 * pi)) * (1 + (u(i)^2 / nu))^(-(nu - 1) / 2);
    rho_2 = (1 / (2 * sigma^2)) * (1 / (2 * pi)^(3/2)) * (gamma((nu + 1) / 2) / (sqrt(nu / 2) * gamma(nu / 2))) * u(i) * (1 + (u(i)^2 / nu))^(-(nu - 1) / 2);

    % Compute the expected Euler characteristic
    EC(i) = mu_0 * rho_0 + mu_1 * rho_1 + mu_2 * rho_2;
end
end
