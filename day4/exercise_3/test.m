clear; clc;

% This sets the precision of the Taylor model
eps_taylor = 1e-6;
% This sets the break condition for the optimization
eps = 1e-6;

% Definition of the given functions
n = 2;
f = @(x) sum(x.^2);
f_start = (1:n)';
g = @(x) sin(x);
g_start = 2;
h = @(x) exp(x(1)/5 + x(2)/2) + x(1)^2 + x(2)^2;
h_start = [5, 3]';
r = @(x) (1-x(1))^2 + 100*(x(2)-x(1)^2)^2;
r_start = [0, 0]';

% Choose the function and starting point
fn = r;
start = r_start;

% Initialize history. Only needed for plotting
x_hist = start';

% Create the initial model
model = firstOrderTaylor(fn, start, eps_taylor);

% Code for plotting the 2D functions
% figure
% fcontour(@(x,y)[fn([x y])], [-0.5 1.5 -0.2 1.5])
% hold on
% pl = plot(x_hist(:,1), x_hist(:,2), '-*', 'Color','blue');
% pl.XDataSource = 'x_hist(:,1)';
% pl.YDataSource = 'x_hist(:,2)';

% Code for plotting 1D
% figure
% fplot(fn, [1, 6])
% hold on 
% pl = plot(x_hist, fn(x_hist), '-*', 'Color','red');
% pl.XDataSource = 'x_hist';
% pl.YDataSource = 'fn(x_hist)';

% Set the maximum number of iterations
k_max = 10000;
% Set the initial step size
delta = 0.1;

k = 0; x = start; found = false;
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
        x = x_new % Update the current point
        % x_hist = [x_hist; x'];
        % refreshdata
        % drawnow
        % pause(0.2)
        model = firstOrderTaylor(fn, x, eps_taylor); % Update the model
    end
    if delta < eps
        found = true;
    end
end

x
% hold on
% refreshdata
% drawnow
% plot(x(1), x(2), 'r*', 'Color','red')
