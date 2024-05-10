clear; clc; close all;

num_cities = 50;
grid_size = 100;

num_generations = 10000;
num_individuals = 100;
num_parents = ceil(num_individuals/10);
p_c = 0.8;
p_m_0 = 0.2;
p_m = p_m_0;

cities = [randi(grid_size, num_cities, 1), randi(grid_size, num_cities, 1)];

P = initializePopulation(cities, num_individuals);
for i = 1:num_individuals
    P{i, 2} = evaluatePath(cities, P{i, 1});
end
P = sortrows(P, 2);
figure
plotPath(gca, cities, P{1, 1}, P{1, 2}, grid_size);

k = 0;
while k < num_generations
    % choose parents
    parents = P(1:num_parents, 1);
    % get all possible pairs of parents
    parents_idx = nchoosek(1:num_parents, 2);
    % create offspring
    for i = 1:height(parents_idx)
        parent1 = parents{parents_idx(i, 1)};
        parent2 = parents{parents_idx(i, 2)};
        O{i, 1} = crossover(parent1, parent2);
        O{i, 2} = evaluatePath(cities, O{i, 1});
    end
    
    P = [P; O];
    % mutation
    for i = 1:num_individuals
        if rand < p_m
            P{i, 1} = randomSwap(P{i, 1});
            P{i, 2} = evaluatePath(cities, P{i, 1});
        end
    end
    P = sortrows(P, 2);
    P = P(1:num_individuals, :);
    
    S = std(cell2mat(P(1:20, 2)));
    if S < 50
        p_m = min(p_m*2, 0.9);
    else
        p_m = max(p_m/2, 0.1);
    end
    
    plotPath(gca, cities, P{1, 1}, P{1, 2}, grid_size);
    drawnow limitrate;
    fprintf('Generation: %d, Best path length: %f, p_m: %f, sigma: %f\n', k, P{1, 2}, p_m, S);
    k = k + 1;
end

function P = initializePopulation(cities, num_individuals)
num_cities = size(cities, 1);
P = cell(num_individuals, 2);
for i = 1:num_individuals
    P{i,1} = randperm(num_cities);
    P{i,2} = evaluatePath(cities, P{i,1});
end
end

function length = evaluatePath(cities, path)
length = 0;
for i = 1:width(path) - 1
    length = length + norm(cities(path(i), :) - cities(path(i + 1), :));
end
% add the distance from the last city back to the first city
% length = length + norm(cities(path(end), :) - cities(path(1), :));
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

function child = crossover(parent1, parent2)
length = width(parent1);
% generate random crossover point
crossover_point = randi([1, length - 1]);
p1slice = parent1(1:crossover_point);
p2slice = parent2(crossover_point+1:end);
child = [p1slice, p2slice];
all = 1:length;
% find missing cities
missing = setdiff(all, child);
% find indices of duplicate cities
[v, w] = unique( child, 'stable' );
duplicate_indices = setdiff( 1:numel(child), w );
% replace duplicate cities with missing cities
child(duplicate_indices) = missing;
if width(unique(child)) ~= length
    error('Crossover failed');
end
end

function child = crossover2(parent1, parent2)
length = width(parent1);
% generate random crossover point
crossover_point = randi([1, length - 1]);
child = [parent1(1:crossover_point+1), parent2];
child = unique(child, 'stable');
if width(child) ~= length
    error('Crossover failed');
end
end

function newPath = randomSwap(oldPath)
% copy the old path to the new path
newPath = oldPath(1:end-1);

% generate two random indices
i = randi([1, size(newPath, 2)]);
j = randi([1, size(newPath, 2)]);

% swap the cities at the two indices
newPath([i, j]) = newPath([j, i]);
newPath = [newPath, newPath(1)];
end

function R = selectParents(P)
N = length(P);
R = cell(N, 1);
for i = 1:2*N
    c = ceil(N*rand(1,2));
    if P{c(1), 2} < P{c(2), 2}
        R{i} = P{c(1), 1};
    else
        R{i} = P{c(2), 1};
    end
end
end

function R = selection(P,F)
N = length(P);
dim = height(P);
R = zeros(dim,2*N);
for j = 1:2*N
    c = ceil(N*rand(1,2));
    if F(c(1)) < F(c(2))
        R(:,j) = P(:,c(1));
    else
        R(:,j) = P(:,c(2));
    end
end
end
