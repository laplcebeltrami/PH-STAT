function PH_boundary_display(X, B)
%   PH_boundary_display(X, B) visualizes the simplicial complex 
%   using only vertex coordinates X and the boundary matrix B. 
%   The function displays:
%       - Vertices as black dots
%       - Directed edges with arrows (gray)
%       - 2-simplices (triangles) in yellow
%       - 3-simplices (tetrahedra) in blue
%
%   Inputs:
%       X - p x d matrix of coordinates of scatter points in R^d
%       B - cell array of boundary matrices defining simplicial relationships
%
%   (C) 2025 Moo K. Chung, 
%   University of Wisconsin-Madison 
%   Email: mkchung@wisc.edu

k = length(B); % Number of dimensions in the simplicial complex

figure;
scatter3(X(:,1), X(:,2), X(:,3), 30, 'k', 'filled'); % Plot vertices
hold on;

% --- Plot directed edges with only arrowheads ---
if k >= 1
    B1 = B{1}; % Boundary matrix for edges
    [v1_idx, edge_idx] = find(B1 == -1); % Source vertices
    [v2_idx, ~] = find(B1 == 1); % Target vertices

    for i = 1:length(edge_idx)
        v1 = X(v1_idx(i), :);
        v2 = X(v2_idx(i), :);
        plot3([v1(1), v2(1)], [v1(2), v2(2)], [v1(3), v2(3)], 'k', 'LineWidth', 1.5); % Draw edge

        % Compute midpoint and direction vector
        mid = (v1 + v2) / 2; % Exact midpoint
        dir_vector = (v2 - v1) / norm(v2 - v1) * 0.05; % Reduced arrowhead length

        % Compute perpendicular vector to form a small triangle
        normal_vector = [0, 0, 1]; % Arbitrary normal vector
        perp_vector = cross(dir_vector, normal_vector);
        if norm(perp_vector) == 0
            perp_vector = cross(dir_vector, [1, 0, 0]); % Use another normal if parallel
        end
        perp_vector = perp_vector / norm(perp_vector) * 0.02; % Reduce width by 1/3

        % Arrowhead triangle vertices
        arrow_tip = mid + dir_vector * 0.5; % Arrow tip slightly forward
        arrow_side1 = mid - perp_vector;
        arrow_side2 = mid + perp_vector;

        % Draw smaller arrowhead as a tiny triangle
        patch([arrow_tip(1), arrow_side1(1), arrow_side2(1)], ...
              [arrow_tip(2), arrow_side1(2), arrow_side2(2)], ...
              [arrow_tip(3), arrow_side1(3), arrow_side2(3)], 'k', 'FaceColor', 'k');
    end
end
% --- Plot 2-simplices (triangles) ---
triangle_vertices_list = [];
if k >= 2
    B2 = B{2}; % Boundary matrix for triangles
    for i = 1:size(B2, 2)
        [edges_idx, ~] = find(B2(:, i)); % Get edges forming a triangle
        vertices = unique([v1_idx(edges_idx); v2_idx(edges_idx)]); % Extract vertices from edge indices
        if length(vertices) == 3
            patch(X(vertices, 1), X(vertices, 2), X(vertices, 3), 'y', 'FaceAlpha', 0.5);
            triangle_vertices_list = [triangle_vertices_list; vertices(:)']; % Store for tetrahedra
        end
    end
end

% --- Plot 3-simplices (tetrahedra) ---
if k >= 3
    B3 = B{3}; % Boundary matrix for tetrahedra
    for i = 1:size(B3, 2)
        [triangles_idx, ~] = find(B3(:, i)); % Get triangles forming a tetrahedron
        tetrahedron_vertices = unique(triangle_vertices_list(triangles_idx, :)); % Extract vertices

        if length(tetrahedron_vertices) == 4
            % Plot tetrahedron faces
            patch(X(tetrahedron_vertices([1 2 3]), 1), X(tetrahedron_vertices([1 2 3]), 2), X(tetrahedron_vertices([1 2 3]), 3), 'b', 'FaceAlpha', 0.3);
            patch(X(tetrahedron_vertices([1 2 4]), 1), X(tetrahedron_vertices([1 2 4]), 2), X(tetrahedron_vertices([1 2 4]), 3), 'b', 'FaceAlpha', 0.3);
            patch(X(tetrahedron_vertices([1 3 4]), 1), X(tetrahedron_vertices([1 3 4]), 2), X(tetrahedron_vertices([1 3 4]), 3), 'b', 'FaceAlpha', 0.3);
            patch(X(tetrahedron_vertices([2 3 4]), 1), X(tetrahedron_vertices([2 3 4]), 2), X(tetrahedron_vertices([2 3 4]), 3), 'b', 'FaceAlpha', 0.3);
        end
    end
end

% Set plot properties
axis equal;
box on;
grid on;
view(3);
set(gca, 'FontSize', 18);
set(gcf, 'Color', 'w', 'InvertHardcopy', 'off');
hold off;
end