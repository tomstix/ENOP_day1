% clear; clc; close all;

% this variable is used to set the number of cities
numCities = 50;
% this variable is used to set the size of the grid
gridSize = 100;

% Algorithm setup
k_max = 1000;
k = 0;
numevals = 0;
p_m = 0.1;
p_c = 0.8;
num_parents = 10;
N = 80;
P = initialize_population(N, numCities, gridSize, cities);
numevals = numevals + N;
P = sortrows(P, 2);

% create the initial plot
path = P{1, 1};
length = P{1, 2};
min_length = length;
figure
p = plot(path(:, 1), path(:, 2), 'o-');
% store gca so we can change the title later
ha = gca;
xlim([0, gridSize])
ylim([0, gridSize])
title(['Path length: ', num2str(length)])
% set the data sources for the plot so we can update the plot later
p.XDataSource = 'path(:, 1)';
p.YDataSource = 'path(:, 2)';
drawnow limitrate

while k < k_max
    
    % crossover
    % select the parents
    parents = P(1:num_parents, 1);
    % combine parents to pairs
    pairs_idx = nchoosek(1:num_parents, 2);
    % select num_parents pairs randomly
    pairs_idx = pairs_idx(randperm(size(pairs_idx, 1), num_parents), :);
    for i = 1:N
        % select parents randomly
        pair = pairs_idx(randi(num_parents), :);
        % select the parents
        parent1 = parents{pair(1)};
        parent2 = parents{pair(2)};
        % crossover
        if rand < p_c
            child = crossover(parent1, parent2);
        else
            child = parent1;
        end
        P{i, 1} = child;
    end
    
    % mutation
    for i = 1:N
        if rand < p_m
            P{i, 1} = randomSwap(P{i, 1});
        end
    end
    
    % evaluate the population
    for i = 1:N
        [c, l] = evaluateGraph(P{i, 1});
        numevals = numevals + 1;
        P{i, 2} = l;
    end
    P = sortrows(P, 2);
    
    % update the best path
    path = [P{1, 1} ; P{1, 1}(1, :)];
    length = P{1, 2};
    if length < min_length
        min_length = length;
        % update the plot
        refreshdata
        drawnow limitrate
        title(['Path length: ', num2str(length)])
    end
    
    fprintf('Generation: %d, Best path length: %f, Function evaluations: %d\n', k, length, numevals);
    
    k = k + 1;
end

function P = initialize_population(N, numCities, gridSize, cities)
P = cell(N, 2);
if nargin < 4
    cities = [randi(gridSize, numCities, 1), randi(gridSize, numCities, 1)];
end
for i = 1:N
    cityOrder = randperm(numCities);
    cityOrder(find(cityOrder == 1, 1)) = [];
    individualcities = zeros(numCities, 2);
    individualcities(1,:) = cities(1,:);
    for j = 2:numCities
        individualcities(j,:) = cities(cityOrder(j - 1),:);
    end
    [individualPath, individualLength] = evaluateGraph(individualcities);
    P(i,:) = {individualPath, individualLength};
end
end

function child = crossover(parent1, parent2)
length = size(parent1, 1);
% generate random crossover point
crossover_point = randi([1, length - 1]);
% create children
child = [parent2(1:crossover_point, :); parent1];
child = unique(child, 'rows', 'stable');
end

function newPath = randomSwap(oldPath)
% copy the old path to the new path
newPath = oldPath;

% generate two random indices that are not the first and the last one
i = randi([2, size(newPath, 1) - 1]);
j = randi([2, size(newPath, 1) - 1]);

% swap the cities at the two indices
newPath([i, j], :) = newPath([j, i], :);
end