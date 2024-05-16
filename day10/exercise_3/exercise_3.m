clear; clc; close all;

lb = [0, 0]; % lower bound for x and y
ub = [5, 5]; % upper bound for x and y

N_0 = 10; % number of initial points
num_iterations = 50; % number of points to add
num_points = N_0 + num_iterations; % total number of points
pause_time = 0.1; % pause time for the plot

x = zeros(num_points, 2);
% add four corners
x(1, :) = lb;
x(2, :) = [ub(1), lb(2)];
x(3, :) = ub;
x(4, :) = [lb(1), ub(2)];
% add random initial points using LHS
x(5:N_0, :) = HypercubeSampling(N_0 - 4, 2) .* (ub - lb) + lb;

% create delaunay triangulation
DT = delaunayTriangulation(x(1:N_0, :));
% plot
figure;
hold on;
points_plot = scatter(x(1:N_0, 1), x(1:N_0, 2), 'filled', 'MarkerFaceColor', '#0072BD');
points_plot.XDataSource = 'x(1:k+1, 1)';
points_plot.YDataSource = 'x(1:k+1, 2)';
tr_plot = triplot(DT,'LineWidth', 0.5, 'Color', '#808080');
hold off

for k = N_0:N_0+num_iterations
    % calculate sizes of the triangles
    tri = DT.ConnectivityList;
    tri_area = zeros(size(tri, 1), 1);
    for i = 1:size(tri, 1)
        tri_area(i) = polyarea(x(tri(i, :), 1), x(tri(i, :), 2));
    end
    % find the biggest triangle and its centroid
    max_area = max(tri_area);
    max_idx = find(tri_area == max_area);
    centroid = incenter(DT, max_idx);
    % plot the new point
    hold on;
    centroid_plot = scatter(centroid(1), centroid(2), 'filled', 'MarkerFaceColor', '#D95319');
    title(sprintf("Points placed: %d of %d", k, num_points))
    fprintf("New Point at: (%f, %f)\n", centroid(1), centroid(2));
    fprintf("Points placed: %d of %d\n\n", k, num_points);
    hold off;
    pause(pause_time);
    % add new point to x
    x(k + 1, :) = centroid;
    % update DT
    DT = delaunayTriangulation(x(1:k + 1, :));
    % update plot
    refreshdata(points_plot, 'caller');
    delete(tr_plot);
    hold on;
    tr_plot = triplot(DT,'LineWidth', 0.5, 'Color', '#808080');
    hold off;
    drawnow;
    pause(pause_time);
    delete(centroid_plot);
end