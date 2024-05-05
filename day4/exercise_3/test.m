clear; clc;

eps_taylor = 1e-6;
eps = 1e-6;

n = 2;
f = @(x) sum(x.^2);
f_start = (1:n)';
g = @(x) sin(x);
g_start = 2;
h = @(x) exp(x(1)/5 + x(2)/2) + x(1)^2 + x(2)^2;
h_start = [5, 3]';
r = @(x) (1-x(1))^2 + 100*(x(2)-x(1)^2)^2;
r_start = [0, 0]';

fn = f;
start = f_start;

model = firstOrderTaylor(fn, start, eps_taylor);

% Code for plotting the 2D functions
% figure
% fcontour(@(x,y)[fn([x y])])

k_max = 10000;
k = 0; x = start; found = false;
delta = 0.1;
while ~found && k <= k_max
    k = k + 1;
    x_new = linearSearch(model, x, delta);
    p = x_new - x;
    actual_reduction = fn(x) - fn(x_new);
    predicted_reduction = fn(x) - model(x_new);
    r = actual_reduction / predicted_reduction;
    
    fprintf('k = %d, r = %d, delta = %d\n', k, r, delta)
    
    if r < 0.25
        delta = 0.25*delta;
    elseif r > 0.75
        delta = 2*delta;
    end
    if r > 0
        x = x_new
        % hold on
        % plot(x(1), x(2), 'r*')
        % hold off
        model = firstOrderTaylor(fn, x, eps_taylor);
    end
    if delta < eps
        found = true;
    end
end

x
% hold on
% plot(x(1), x(2), 'r*', 'Color','blue')
