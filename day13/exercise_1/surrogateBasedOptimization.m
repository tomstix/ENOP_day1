clear; clc; close all;

% Function definition
func = @(x) exercise_1_function(x);
lb = [-2, -1]; % lower bound for x and y
ub = [3, 2]; % upper bound for x and y

% Setup variables
N_0 = 25; % number of initial points
% Create uniform grid
points_per_dim = sqrt(N_0);
if mod(points_per_dim, 1) ~= 0
    error('N must be a perfect square');
end
xs = linspace(lb(1), ub(1), points_per_dim)';
ys = linspace(lb(2), ub(2), points_per_dim)';

% Create/load sample points from the funciton
if exist("samples.mat", 'file')
    load("samples.mat");
else
    % Create all combinations of grid points
    [p, q] = meshgrid(xs, ys);
    combs = [p(:) q(:)];

    % Evaluate function at grid points
    zs = arrayfun(@(i) func(combs(i,:)), 1:size(combs, 1));
    
    % Save samples
    save("samples.mat", "combs", "zs");
end

% create Phi matrix
Phi = zeros(N, N);
for i = 1:N
    for j = 1:N
        Phi(i, j) = phi_r(norm(combs(:,i) - combs(:,j)));
    end
end
    
% solve linear system
lambda = Phi \ zs';
    
% create surrogate model
sunc = @(x) sum(arrayfun(@(i) lambda(i).*phi_r(norm(x - combs(:,i))), 1:N));



% % Low-fidelity model created with datafitting
% parameters = ones(1, N_0 + 1);
% cunc = @(a, x) sum(a(1)*x.^10 + a(2)*x.^9 + a(3)*x.^8 + a(4)*x.^7 + a(5)*x.^6 + a(6)*x.^5 + a(7)*x.^4 + a(8)*x.^3 + a(9)*x.^2 + a(10)*x + a(11));
% 
% % Find function that takes a, b, c and d as parameters and returns the
% % difference between the model and the actual data
% d = samplePoints;
% y = sampleValues;
% fun = @(x) cunc(x, d) - y;
% 
% % Optimize parameters
% parameters = lsqnonlin(fun, parameters);
% 
% % Update low-fidelity model
% cunc = @(x) cunc(parameters,x);
% 
% fcontour((@(x,y) cunc([x y])));
% 
% % Optimization
% x_0 = [0 0];
% 
% % Tangent function and gradient at x^(0) using finite differences.
% [func_, func_grad_x_0] = firstOrderTaylor(func, x_0, 1e-12);
% [cunc_, cunc_grad_x_0 ] = firstOrderTaylor(cunc, x_0, 1e-12);
% % func_grad_x_0 = func_(x_0);
% % cunc_grad_x_0 = cunc_(x_0);
% 
% % Beta-correlation method
% funcx_0 = func(x_0);
% cuncx_0 = cunc(x_0);
% aunc = @(x) ((funcx_0 / cuncx_0) + ((func_grad_x_0 * cuncx_0 - cunc_grad_x_0 * funcx_0)/cuncx_0^2) .* (x - x_0));
% sunc = @(x) aunc(x).*cunc(x);
% 
% fcontour((@(x,y) sunc([x y])));