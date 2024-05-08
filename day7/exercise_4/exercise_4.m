clear; clc;

f = @auckley; % set the function to optimize
dim = 2; % dimensionality of the search space
lim = [-2, 2]; % parameter space limits

max_iterations = 300; % maximum number of iterations
num_particles = 100; % number of particles
chi = 0.7298; % inertia weight
c1 = 2.05; % cognitive weight
c2 = 2.05; % social weight

drawnow limitrate
[f_best, g] = particle_swarm_optimization(f, dim, num_particles, lim, max_iterations, chi, c1, c2, true);