% Function to plot the path on the grid. This version only uses the path
function plotPath_alt(ax, path, length, grid_size)
cla(ax);
plot(ax, path(:, 1), path(:, 2), 'o');
hold(ax, 'on');
plot(ax, path(:, 1), path(:, 2), 'r-');
xlim(ax, [0, grid_size]);
ylim(ax, [0, grid_size]);
title(ax, ['Path length: ', num2str(length)]);
end