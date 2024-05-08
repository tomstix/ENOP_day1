
clear; clc;

max_iterations = 1000; % maximum number of iterations
lambda_0 = 0.5; % lambda_0 parameter for smart random search
num_runs = 100; % number of runs for the experiment

% uncomment the following lines to run a single optimization
fn = @rosenbrock; % function to optimize
lim = [-2 2]; % parameter space limits
dim = 2; % dimensionality of the search space
f_best = random_search(fn, lim, dim, max_iterations, true)
f_best_s = smart_random_search(fn, lim, dim, max_iterations, lambda_0, true)

% or run an experiment to compare the two alogirthms with all the functions
% run_experiment(num_runs, max_iterations, lambda_0);

function run_experiment(num_runs, max_iterations, lambda_0)
for function_select = 1:3
    switch function_select
        case 1
            fn = @rosenbrock;
            lim = [-2 2];
            dims = 2;
            name = 'Rosenbrock';
        case 2
            fn = @fp;
            lim = [-2 2];
            dims = [1,2];
            name = 'f_p';
        case 3
            fn = @auckley;
            lim = [-10 10];
            dims = [1,2,3];
            name = 'Auckley';
    end
    for j = 1:width(dims)
        dimensionality = dims(j);
        f_best = zeros(num_runs, 1);
        f_best_s = zeros(num_runs, 1);
        for i = 1:num_runs
            f_best(i) = random_search(fn, lim, dimensionality, max_iterations, false);
            f_best_s(i) = smart_random_search(fn, lim, dimensionality, max_iterations, lambda_0, false);
        end
        disp(['Function: ', name, ' (n = ', num2str(dimensionality), ')'])
        disp(['Random search - Best: ', num2str(min(f_best)), ' - Worst: ', num2str(max(f_best)), ' - Mean: ', num2str(mean(f_best)), ' - Std: ', num2str(std(f_best))])
        disp(['Smart random search - Best: ', num2str(min(f_best_s)), ' - Worst: ', num2str(max(f_best_s)), ' - Mean: ', num2str(mean(f_best_s)), ' - Std: ', num2str(std(f_best_s))])
        fprintf('\n')
    end
end
end

function x = generate_random_point(lim, dimensionality)
% generates a random vector within the limits specified by lim and size specified by dimensionality
x = lim(1) + (lim(2) - lim(1)) * rand(dimensionality, 1);
end

function [f_best, x_best] = random_search(fn, lim, dimensionality, k_max, print)
f_best = inf; % initialize the best value found so far
k = 0; % initialize the iteration counter
num_evaluations = 0; % initialize the number of function evaluations
while k < k_max
    k = k + 1;
    % choose a random point in the search space
    x = generate_random_point(lim, dimensionality);
    % evaluate the function at that point
    f = fn(x);
    num_evaluations = num_evaluations + 1;
    % update the best solution found so far
    if f < f_best
        f_best = f;
        x_best = x;
    end
    if print
        disp(['Iteration: ', num2str(k), ' - Best value: ', num2str(f_best), ' - Evaluations: ', num2str(num_evaluations)])
    end
end
end

function [f_best, x_best] = smart_random_search(fn, lim, dimensionality, k_max, lambda_0, print)
x_best = generate_random_point(lim, dimensionality); % intial random point
f_best = fn(x_best); % initialize best value with the value at the initial point
k = 0;
update_lambda = @(k, k_max, lambda_0) lambda_0 * (1 - k/k_max); % update rule for lambda
lambda = lambda_0; % initialize lambda
num_evaluations = 0;
while k < k_max
    k = k + 1;
    % choose a random point in the search space
    x = generate_random_point(lim, dimensionality);
    x = lambda * x + (1 - lambda) * x_best;
    % evaluate the function at that point
    f = fn(x);
    num_evaluations = num_evaluations + 1;
    % update the best solution found so far
    if f < f_best
        f_best = f;
        x_best = x;
    end
    lambda = update_lambda(k, k_max, lambda_0);
    if print
        disp(['Iteration: ', num2str(k), ' - Best value: ', num2str(f_best), ' - Evaluations: ', num2str(num_evaluations), ' - Lambda: ', num2str(lambda)])
    end
end
end