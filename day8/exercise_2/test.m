% This script is used to compare the two evolutionary algorithms with the random search algorithm.
clear; clc; close all;

num_cities = randi([20 500]); % Number of cities
grid_size = 100; % Size of the grid

num_individuals = 50; % Number of individuals for the evolutionary algorithms
max_evals = 500000; % computational budget
p_c = 0.8; % Crossover probability for the evolutionary algorithms
p_m = 1/num_individuals; % Mutation probability for the evolutionary algorithms

% Generate random cities
cities = [randi(grid_size, num_cities, 1), randi(grid_size, num_cities, 1)];

figure
[~, initial_length] = evaluateGraph(cities);
plotPath(gca, cities, 1:num_cities, initial_length, grid_size);

% Select which algorithm to run
% evolutionaryTSP(cities, max_evals, num_individuals, p_c, p_m, true, true, grid_size);
evolutionaryTSP2(cities, num_individuals, max_evals, p_c, p_m, true, true, grid_size);
% randomTSP(cities, max_evals, true, true, grid_size)
