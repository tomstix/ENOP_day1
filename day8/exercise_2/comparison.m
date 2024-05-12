% This script is used to compare the two evolutionary algorithms with the random search algorithm.
clear; clc; close all;

num_cities = 20; % Number of cities
grid_size = 100; % Size of the grid

num_individuals = 100; % Number of individuals for the evolutionary algorithms
max_evals = 500000; % computational budget
p_c = 0.8; % Crossover probability for the evolutionary algorithms
p_m = 0.2; % Mutation probability for the evolutionary algorithms

num_runs = 16; % Number of runs to average over

% Generate random cities that are used for all runs
cities = [randi(grid_size, num_cities, 1), randi(grid_size, num_cities, 1)];

% Run the algorithms num_runs times (in parallel to speed up the process)
parfor run = 1:num_runs
    fprintf('Run %d\n', run)
    history1{run} = evolutionaryTSP(cities, max_evals, num_individuals, p_c, p_m);
    history2{run} = evolutionaryTSP2(cities, num_individuals, max_evals, p_c, p_m, false, false, grid_size);
    historyR{run} = randomTSP(cities, max_evals, false, false, grid_size);
end

% Interpolate and average the results
xs = 1:max_evals;
parfor run = 1:num_runs
    ys1(:, run) = interp1(history1{run}(:, 2), history1{run}(:, 3), xs, 'previous');
    ys2(:, run) = interp1(history2{run}(:, 2), history2{run}(:, 3), xs, 'previous');
    ysR(:, run) = historyR{run}(:, 3);
end
history1 = [xs', mean(ys1, 2)];
history2 = [xs', mean(ys2, 2)];
historyR = [xs', mean(ysR, 2)];
history1 = rmmissing(history1);
history2 = rmmissing(history2);

% Plot the results
figure
hold on
plot(history1(:, 1), history1(:, 2), 'DisplayName', 'evolutionaryTSP')
plot(history2(:, 1), history2(:, 2), 'DisplayName', 'evolutionaryTSP2')
plot(historyR(:, 1), historyR(:, 2), 'DisplayName', 'randomTSP')
title(sprintf("Max evals: %.1e, num individuals: %d, p_c: %.2f, p_m: %.2f", max_evals, num_individuals, p_c, p_m))
xlabel('Number of evaluations')
ylabel('Best length')
xscale log
xlim([num_individuals, max_evals])
legend