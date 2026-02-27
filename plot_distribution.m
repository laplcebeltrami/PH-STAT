function plot_distribution(stat_s, bins, observation, binRange, bincolor)
% function plot_distribution(stat_s, bins, observation, binRange, bincolor)
%
% INPUTS:
%   stat_s:      data vector (e.g., null distribution)
%   bins:        number of histogram bins
%   observation: observed value(s) to be plotted as vertical line(s)
%   binRange:    [optional] 2-element vector specifying min and max bin range
%   bincolor:    [optional] 3-element RGB vector for histogram and line color
%
% The code is part of the PH-STAT (Statistical Inference on Persistent Homology) package
% and is available at https://github.com/laplcebeltrami/PH-STAT
%
% (C) 2022 D. Vijay Anand, Moo K. Chung (mkchung@wisc.edu)
% University of Wisconsin-Madison
%
% Update history:
%   2022 - Created (Anand & Chung)
%   2023 - Added optional arguments for binRange and bincolor to handle variable
%          number of input arguments gracefully.
%   2025 argument based error fixed
%   2026 Feb 26, binrange fixed to account for extrem rate events. 

% If fewer than 5 arguments are provided, define defaults for binRange and bincolor.
if nargin < 5 || isempty(bincolor)
    bincolor = [0.7, 0.7, 0.7];  % default gray
end

if nargin < 4 || isempty(binRange)
    % Robust default range: avoid being dominated by extreme outliers.
    % Use central quantiles (winsorized range) so histogram does not collapse. 
    stat_s0 = stat_s(isfinite(stat_s));
    if isrow(stat_s0), stat_s0 = stat_s0'; end

    if isempty(stat_s0)
        binRange = []; % no data  
    else
        % central 0.5%–99.5% range (adjust if you want)
        qL = prctile(stat_s0, 0.5);
        qU = prctile(stat_s0, 99.5);

        if qU == qL
            % near-constant distribution: create a small window
            eps0 = 1e-6 + abs(qL)*1e-6;
            binRange = [qL - eps0, qU + eps0];  
        else
            % small padding so edge values are not stuck on borders
            pad = 0.02*(qU - qL);
            binRange = [qL - pad, qU + pad];    
        end
    end
end


% Ensure stat_s is a column vector
[rows, cols] = size(stat_s);
if rows == 1
    stat_s = stat_s';
end

% Plot histogram
if ~isempty(binRange)
    binEdges = linspace(binRange(1), binRange(2), bins + 1);
    histogram(stat_s, binEdges, 'FaceColor', bincolor, 'Normalization', 'probability');
else
    histogram(stat_s, bins, 'FaceColor', bincolor, 'Normalization', 'probability');
end

% Optionally plot vertical line(s) for observed values
if ~isempty(observation)
    hold on;
    % observation might be scalar or vector; plot each element
    for obsVal = observation
        plot([obsVal, obsVal], [0, 1], '--', 'LineWidth', 2, 'Color', 'r');
    end
end

% Aesthetics and figure properties
set(gcf, 'Position', [400 400 600 250]);
set(gca, 'FontSize', 16);
xlabel('Test Statistic');
ylabel('Probability');
title('Distribution');


