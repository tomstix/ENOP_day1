function visualizePath(ax, points, isFinal)
% Plot the path
plot(ax, points(:,1), points(:,2), 'black', 'LineWidth', 1);
if height(points) > 1
    plot(ax, points(1,1), points(1,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    if isFinal
        plot(ax, points(end,1), points(end,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    end
end
drawnow
end