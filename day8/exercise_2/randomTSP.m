function [history, final_P] = randomTSP(cities, max_evals, print, plot, grid_size)
if nargin < 3
    print = false;
end
if nargin < 4
    plot = false;
end
num_cities = size(cities, 1);
perm = randperm(num_cities);
best_path = perm;
best_length = evaluatePath(cities, best_path);
history = [];
if plot
    figure
    plotPath(gca, cities, best_path, best_length, grid_size);
    drawnow
end
k = 0;
while k < max_evals
    k = k + 1;
    perm = randomSwap(best_path);
    length = evaluatePath(cities, perm);
    if length < best_length
        best_length = length;
        best_path = perm;
        history = [history; k, k, best_length];
    end
    if plot
        plotPath(gca, cities, best_path, best_length, grid_size);
        drawnow limitrate
    end
    if print
        fprintf('Iteration: %d, Best length: %f\n', k, best_length);
    end
end
history = [history; k, k, best_length];
final_P = cities(best_path, :);
end

function length = evaluatePath(cities, path)
length = 0;
for i = 1:width(path) - 1
    length = length + norm(cities(path(i), :) - cities(path(i + 1), :));
end
% add the distance from the last city back to the first city
length = length + norm(cities(path(end), :) - cities(path(1), :));
end

function newPath = randomSwap(oldPath)
% copy the old path to the new path
newPath = oldPath;

% generate two random indices
i = randi([1, size(newPath, 2)]);
j = randi([1, size(newPath, 2)]);

% swap the cities at the two indices
newPath([i, j]) = newPath([j, i]);
end

function plotPath(ax, cities, path, length, grid_size)
cla(ax);
plot(ax, cities(:, 1), cities(:, 2), 'o');
hold(ax, 'on');
plot(ax, cities([path, path(1)], 1), cities([path, path(1)], 2), 'r-');
xlim(ax, [0, grid_size]);
ylim(ax, [0, grid_size]);
title(ax, ['Path length: ', num2str(length)]);
end