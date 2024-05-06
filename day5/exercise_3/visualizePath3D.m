function visualizePath3D(ax, points, fn, isFinal)
points(:,3) = fn(points(:,1), points(:,2)); % calculate the z values
% Plot the path
plot3(ax, points(:,1), points(:,2), points(:,3), 'r', 'LineWidth', 1);
if height(points) > 1
    % plot the start point
    plot3(ax, points(1,1), points(1,2), points(1,3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    if isFinal
        % plot the end point
        plot3(ax, points(end,1), points(end,2), points(end,3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    end
end
drawnow
end