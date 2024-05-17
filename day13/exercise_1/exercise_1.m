clear; clc; close all;

addpath('./dace/');

dim = 2; % dimensionality of the search space
lim = [-2, 3; -1, 2]; % parameter space limits
fn = @exercise_1_function;
k_max = 30; % maximum number of iterations
eps = 1e-6; % convergence criterion

if exist('samples.mat', 'file') == 2
    load('samples.mat'); % this file contains some presampled points
else
    grid_size = 4;

    % Generate initial samples on uniform grid
    init_samples = gridsamp(lim', [grid_size, grid_size]);
    
    % Evaluate the function
    init_values = zeros(height(init_samples), 1);
    for i = 1:height(init_samples)
        fprintf('Evaluating sample %d/%d\n', i, height(init_samples));
        init_values(i) = fn(init_samples(i,:));
    end
    
    save('samples.mat', 'init_samples', 'init_values');
end

% create the initial model
dmodel = dacefit(init_samples, init_values, @regpoly0, @corrgauss, [1, 1], [0.01, 0.01], [10,10]);

% plot the initial model
xs = linspace(lim(1,1), lim(1,2), 100);
ys = linspace(lim(2,1), lim(2,2), 100);
[xi, yi] = meshgrid(xs, ys);
zi = predictor([xi(:), yi(:)], dmodel);

figure;
ax = gca;
hold on;
sc_samp = scatter(init_samples(:,1), init_samples(:,2), 'r');
cont = contour(xi, yi, reshape(zi, size(xi)), 20);

samples = init_samples;
values = init_values;
last_val = Inf;
f_best = min(values);
x_best = samples(values == f_best, :);
num_nonimprovements = 0; % number of iterations without improvement
for k = 1:k_max
    % optimize the model using particle swarm optimization
    s = @(x) predictor(x, dmodel);
    x0 = [1, 1];
    lb = [lim(1,1), lim(2,1)];
    ub = [lim(1,2), lim(2,2)];
    fprintf('Optimizing model...\n');
    [fval, x] = particle_swarm_optimization(s, dim, 100, lb, ub);
    fprintf('Found minimum of surrogate model at [%f, %f]\n', x(1), x(2));
    
    % plot the minimum point
    sc_min = scatter(x(1), x(2), 'g', 'filled');

    % add the new point to the samples
    fprintf('Evaluating function at [%f, %f]\n', x(1), x(2));
    f_new = fn(x);
    delta = last_val - f_new;
    last_val = f_new;
    if f_new < f_best
        f_best = f_new;
        x_best = x;
        num_nonimprovements = 0;
    else
        num_nonimprovements = num_nonimprovements + 1;
        if num_nonimprovements > 5
            fprintf('Optimization terminated because there was no improvement. Total function evaluations: %d\n', height(values)); 
            break;
        end
    end
    samples = [samples; x];
    values = [values; f_new];

    fprintf('Iteration %d: x = [%f, %f], f(x) = %f, delta = %f\n', k, x(1), x(2), f_new, delta);
    fprintf('Best value: %f at [%f, %f]\n', f_best, x_best(1), x_best(2));
    % break if the change in the function value is less than eps
    if abs(delta) < eps
        fprintf('Optimization terminated because of convergence. Total function evaluations: %d\n', height(values)); 
        break;
    end

    % update the model
    fprintf('Updating model...\n');
    dmodel = dacefit(samples, values, @regpoly0, @corrgauss, [1, 1], [0.01, 0.01], [10,10]);

    % plot the updated model
    zi = predictor([xi(:), yi(:)], dmodel);
    cla(ax);
    cont = contour(xi, yi, reshape(zi, size(xi)), 20);
    sc_samp = scatter(samples(:,1), samples(:,2), 'r');
    drawnow;
end

scatter(x_best(1), x_best(2), 'b', 'filled');
title(sprintf('Initial samples: %d, Total samples: %d, Minimum: %0.6f', height(init_samples), height(samples), f_best));
