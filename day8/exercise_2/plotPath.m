% Function to plot the path. This version uses the city coordinates and the path as a permutation of the city indices
function plotPath(ax, cities, path, length, grid_size)
cla(ax);
plot(ax, cities(:, 1), cities(:, 2), 'o');
hold(ax, 'on');
plot(ax, cities([path, path(1)], 1), cities([path, path(1)], 2), 'r-');
xlim(ax, [0, grid_size]);
ylim(ax, [0, grid_size]);
title(ax, ['Path length: ', num2str(length)]);
end