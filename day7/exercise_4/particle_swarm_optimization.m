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

function [f_best, g, g_hist, f_hist, k, f_calls] = particle_swarm_optimization(f, dim, n, lim, k_max, chi, c1, c2, print)

% function to generate random vectors
generate_random_vectors = @(lim, dim, n) lim(1) + (lim(2) - lim(1)) * rand(dim, n);

if dim == 2 && print % plot the function if it is 2D
    fcontour(@(x,y)f([x;y]), lim);
    hold on
    xlim(lim);
    ylim(lim);
end
% initialize x and v
x = generate_random_vectors(lim, dim, n);
v_up = (lim(2) - lim(1));
v = generate_random_vectors([-v_up, v_up], dim, n);
x_best = x;
fs = f(x);
[f_best, g_idx] = min(fs);
g = x(:,g_idx);
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
        r1 = rand(dim, 1);
        r2 = rand(dim, 1);
        
        % update velocity
        v(:,j) = chi * (v(:,j) + c1 * r1 .* (x_best(:,j) - x(:,j)) + c2 * r2 .* (g - x(:,j)));
        
        % update position
        x(:,j) = x(:,j) + v(:,j);
        
        f_x = f(x(:,j));
        % update local best
        if f_x < fs(j)
            x_best(:,j) = x(:,j);
            fs(j) = f_x;
        end
        if f_x < f_best
            g = x(:,j);
            f_best = f_x;
        end
        f_calls = f_calls + 1; % we had to evaluate the function once
    end
    g_hist(:,k) = g;
    f_hist(k) = f_best;
    if print
        fprintf('Iteration: %d - Best value: %0.4e - Evaluations: %d\n', k, f_best, f_calls)
        if dim == 2
            plot(x_best(1,:), x_best(2,:), 'b*', 'MarkerSize', 0.5)
            plot(g(1), g(2), 'r*', 'MarkerSize', 5)
            drawnow
        end
    end
end
if print
    gst = strjoin(cellstr(num2str(g(:))),', ');
    disp(['Minimum value found: ', num2str(f_best), ' at: [', gst, ']'])
end
if dim == 2 && print
    plot(g(1), g(2), 'r*', 'MarkerSize', 20)
    hold off
end
end