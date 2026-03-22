function graph_overlay_timeseries(Nodes, Z)
% GRAPH_OVERLAY_TIMESERIES overlay normalized time series at node positions
%
%   graph_overlay_timeseries(Nodes, Z)
%
% Inputs
%   Nodes : n x 2 coordinates
%   Z     : n x T matrix of time series (n nodes, T samples)
%
% (c) 2025 Moo K. Chung
% University of Wisconsin-Madison
% mkchung@wisc.edu 
% 
% The code is downloaded from 
% https://github.com/laplcebeltrami/hodge
% If you are using the code, refernce one of Hodge papers listed in GitHub.  
%
% Update history: November 7, 2025


[n, T] = size(Z);
Tshow  = min(200, T);           % show up to 200 timepoints
xgrid  = linspace(-1, 1, Tshow);

% scaling and displacement
xspan  = 0.30;     % width of mini-traces
yscale = 0.15;     % height scale
pushx  = 0.5;      % outward displacement in x
pushy  = 0.2;      % outward displacement in y

% center of graph
ctr = mean(Nodes,1);

hold on;
for i = 1:n
    % normalize series
    ts = Z(i,1:Tshow);
    ts = ts - mean(ts);
    m  = max(abs(ts));
    if m>0, ts = ts/m; end

    % vector from center → node
    v = Nodes(i,:) - ctr;
    if norm(v) > 0
        v = v / norm(v);
    end
    
    % apply separate x/y push
    offset = [pushx*v(1), pushy*v(2)];

    % plot mini time series
    %plot(Nodes(i,1)+offset(1) + xspan*xgrid, ...
    %     Nodes(i,2)+offset(2) + yscale*ts, ...
    %     'Color',[1 0.4 0.4], 'LineWidth',1.2);

    plot(Nodes(i,1)+offset(1) + xspan*xgrid, ...
         Nodes(i,2)+offset(2) + yscale*ts, ...
         'Color','g', 'LineWidth',1);

end
axis equal;
end

