% Perform particle swarm optimization
% f: function to optimize
% dim: dimensionality of the search space
% n: number of particles
% lim: parameter space limits
% k_max: maximum number of iterations
% chi: inertia weight
% c1: cognitive weight
% c2: social weight
% print: whether to print (and plot for 2D) the results

function [f_best, g, g_hist, f_hist, k, f_calls] = particle_swarm_optimization(f, dim, n, lb, ub, k_max, chi, c1, c2, print)
if nargin < 10
    print = false;
end
if nargin < 9
    chi = 0.7298;
end
if nargin < 8
    c1 = 2.05;
end
if nargin < 7
    c2 = 2.05;
end
if nargin < 6
    k_max = 300;
end

% function to generate random vectors
% generate_random_vectors = @(lim, dim, n) lim(1) + (lim(2) - lim(1)) * rand(dim, n);
generate_random_vectors = @(n, dim, lb, ub) lb + ((ub - lb) .* rand(n, dim));

% initialize x and v
x = generate_random_vectors(n, dim, lb, ub);
v_up = (ub - lb);
v = generate_random_vectors(n, dim, -v_up, v_up);
x_best = x;
fs = f(x);
[f_best, g_idx] = min(fs);
g = x(g_idx, :);
g_hist = zeros(dim, k_max);
f_hist = zeros(1, k_max);
g_hist(:,1) = g;
f_hist(1) = f_best;
k = 0; % initialize the iteration counter
f_calls = n; % initialize the number of function evaluations (we have evaluated the function for each particle)
brk = false;
while k < k_max && ~brk
    k = k+1;
    for j = 1:n
        % generate random vectors
        r1 = rand(1, dim);
        r2 = rand(1, dim);
        
        % update velocity
        v(j, :) = chi * (v(j, :) + c1 * r1 .* (x_best(j, :) - x(j, :)) + c2 * r2 .* (g - x(j, :)));
        
        % update position
        x(j, :) = x(j, :) + v(j, :);

        % constrain the position
        x(j, :) = max(x(j, :), lb);
        x(j, :) = min(x(j, :), ub);
        
        f_x = f(x(j, :));
        % update local best
        if f_x < fs(j)
            x_best(j, :) = x(j, :);
            fs(j) = f_x;
        end
        if f_x < f_best
            g = x(j, :);
            f_best = f_x;
        end
        f_calls = f_calls + 1; % we had to evaluate the function once
    end
    g_hist(:,k) = g;
    f_hist(k) = f_best;
    if print
        fprintf('Iteration: %d - Best value: %0.4e - Evaluations: %d\n', k, f_best, f_calls)
    end
end
if print
    gst = strjoin(cellstr(num2str(g(:))),', ');
    disp(['Minimum value found: ', num2str(f_best), ' at: [', gst, ']'])
end
end