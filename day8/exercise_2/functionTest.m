clear
clc
close all

% this variable is used to set the number of cities
numCities = 50;
% this variable is used to set the size of the grid
gridSize = 100;

% create random cities (coordinates) within the grid
cities = [randi(gridSize, numCities, 1), randi(gridSize, numCities, 1)];

% Function parameters:
N = 100; % Population size
Pc = 0.8; % Crossover probability
Pm = 1 / N; % Mutation probability
maxEvaluations = 100000; % Max evaluations

[history, bestPath] = evolutionaryTSP(cities, maxEvaluations, N, Pc, Pm);