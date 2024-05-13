clear; clc;

f = @auckley; % set the function to optimize
dim = 3; % dimensionality of the search space
lim = [-10, 10]; % parameter space limits

num_individuals = 100; % number of individuals
num_generations = 1000; % number of generations
pc = 0.8; % crossover probability
pm0 = 0.1; % initial mutation probability
beta = 3; % mutation parameter

if dim == 2 % plot the function if it is 2D
    fcontour(@(x,y)f([x;y]), lim);
    hold on
    xlim(lim);
    ylim(lim);
end

fga(f, dim, lim, num_individuals, num_generations, pc, pm0, beta, @display_results);

function display_results(k, P, F, pm, evals)
f_best = min(F);
disp(['Generation: ', num2str(k), ', Best fitness: ', num2str(f_best), ', Mutation probability: ', num2str(pm), ', Evaluations: ', num2str(evals)]);

persistent sc;
delete(sc);
if height(P) ~= 2
    return;
end
sc = scatter(P(1,:),P(2,:),50,'r.');
drawnow
end
