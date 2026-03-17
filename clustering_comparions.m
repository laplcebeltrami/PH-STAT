function perf = clustering_comparions(nSim, sigma, npoints, nSubjects, mode)
% clustering_comparions
%
% Run clustering simulations for either equal or different topology.
%
% INPUT
%   nSim    : number of simulations
%   sigma   : noise level
%   npoints : number of scatter points
%   nSubjects: numberf of subjects in each group
%   mode    : 'equal' or 'different'

%
% OUTPUT
%   perf : structure containing mean and std of accuracies
%
% (C) 2026 Moo K. Chung
% Univeristy of Wisconsin-Madison


% Initialize accuracy storage
acc_K  = zeros(nSim,1);
acc_B0 = zeros(nSim,1);
acc_B1 = zeros(nSim,1);
acc_H  = zeros(nSim,1);
acc_WS = zeros(nSim,1);

for i = 1:nSim
    
    % Generate simulated data
    if strcmp(mode,'equal')
        coord = simulate_circle_equal(sigma, npoints, nSubjects);
    else
        coord = simulate_circle_difference(sigma, npoints,nSubjects);
    end
    
    g = coord2dist(coord);

    % k-means clustering
    acc = kmeans_cluster(g);
    acc_K(i) = acc.mean;

    % bottleneck distance
    acc = WS_bottleneck_cluster(g);
    acc_B0(i) = acc.zero;
    acc_B1(i) = acc.one;

    % hierarchical clustering
    acc_H(i) = GH_cluster(g);

    % Wasserstein distance
    acc = WS_cluster(g);
    acc_WS(i) = acc.mean;
end

% Performance metrics
perf.K.mean  = mean(acc_K);  perf.K.std  = std(acc_K);
perf.B0.mean = mean(acc_B0); perf.B0.std = std(acc_B0);
perf.B1.mean = mean(acc_B1); perf.B1.std = std(acc_B1);
perf.H.mean  = mean(acc_H);  perf.H.std  = std(acc_H);
perf.WS.mean = mean(acc_WS); perf.WS.std = std(acc_WS);


% Display performance summary
disp('Performance summary (mean ± std)')
fprintf('K-means           : %.3f ± %.3f\n', perf.K.mean,  perf.K.std);
fprintf('Bottleneck-0      : %.3f ± %.3f\n', perf.B0.mean, perf.B0.std);
fprintf('Bottleneck-1      : %.3f ± %.3f\n', perf.B1.mean, perf.B1.std);
fprintf('Gromov-Hausdorff  : %.3f ± %.3f\n', perf.H.mean,  perf.H.std);
fprintf('Wasserstein       : %.3f ± %.3f\n', perf.WS.mean, perf.WS.std);

end