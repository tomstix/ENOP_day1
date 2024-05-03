clear; clc;

eps = 1e-6;

g = @(x) sin(x);
g_start = 2;
h = @(x) exp(x(:,1)/5 + x(:,2)/2) + x(:,1).^2 + x(:,2).^2;
h_start = [5, 3]';
r = @(x) (1-x(1))^2 + 100*(x(2)-x(1)^2)^2;
r_start = [0, 0]';

fn = h;
start = h_start;
model_x0 = start;
model = firstOrderTaylor(fn, model_x0, eps);

figure
hold on
% fplot(fn, [0, 5], 'b');
% fplot(model, [0, 5], 'r');
fcontour(@(x,y)fn([x, y]));
hold off
ax = gca;

k_max = 100;
k = 0; x = start; delta_min = 1e-6; delta_max = 1; found = false;
delta = delta_max;
while ~found && k <= k_max
    k = k + 1;
    cla(ax);
    hold on
    plot(x);
    hold off
    drawnow
    h_tr = goldenRatioSearch(model, x-delta, x+delta, 1e-6);
    p = h_tr - x;
    fprintf('p = %f\n', p)
    actual_reduction = fn(x) - fn(x + p);
    predicted_reduction = fn(x) - model(x + p);
    fprintf('Actual reduction: %f, Predicted reduction: %f\n', actual_reduction, predicted_reduction);
    if predicted_reduction == 0
        r = 0;
    else
        r = actual_reduction / predicted_reduction;
    end
    
    fprintf('k = %d, x = %f, h_tr = %f, r = %f, delta = %f\n', k, x, h_tr, r, delta)
    
    if r < 0.25
        delta = max(0.25*delta, delta_min);
    elseif r > 0.75 && abs(h_tr) == delta
        delta = min(2*delta, delta_max);
    end
    if r > 0
        x = h_tr;
        model_x0 = x;
        model = firstOrderTaylor(fn, model_x0, eps);
    end
    pause(0.5)
end
