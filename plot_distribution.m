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

% If fewer than 5 arguments are provided, define defaults for binRange and bincolor.
if nargin < 5 || isempty(bincolor)
    bincolor = [0.7, 0.7, 0.7];  % default gray
end
if nargin < 4 || isempty(binRange)
    binRange = [];  % no specific bin range
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

%whitebg(gcf,'w');
set(gcf,'Color','w','InvertHardcopy','off');

