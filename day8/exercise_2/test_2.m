clear; clc; close all;

num_cities = 50;
grid_size = 100;

num_generations = 5000;
num_individuals = 50;
max_evals = 100000;
p_c = 0.8;
p_m_0 = 0.1;
p_m = p_m_0;

cities = [randi(grid_size, num_cities, 1), randi(grid_size, num_cities, 1)];

[history, final_P] = evolutionaryTSP2(cities, num_individuals, num_generations, max_evals, p_c, p_m_0, false, false, grid_size);

[hist2, final_P2] = randomTSP(cities, max_evals, false, false, grid_size);

figure
plot(history(:, 2), history(:, 3), 'DisplayName', 'Evolutionary TSP')
hold on
plot(hist2(:, 2), hist2(:, 3), 'DisplayName', 'Random TSP')
legend
title('Comparison of Evolutionary TSP and Random TSP')
xlabel('Function Evaluations')
ylabel('Path length')