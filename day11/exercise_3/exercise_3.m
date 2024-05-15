clear; clc; close all;

% Gaussian basis function
c = 1;
phi_r = @(r) exp(-c*r.^2);

% N's to test
Ns = [9 25 49 100 169];

% function limits
lim = [-3 3];

% plot true function
figure
fmesh(@(x,y) exercise_3_function([x;y]), lim);
xlim([-4 4])
ylim([-4 4])
title('True function');
% saveas(gcf, 'd11e3_true.png');
drawnow

for N = Ns
    points_per_dim = sqrt(N);
    if mod(points_per_dim, 1) ~= 0
        error('N must be a perfect square');
    end
    
    % create uniform grid
    xs = linspace(lim(1), lim(2), points_per_dim);
    ys = linspace(lim(1), lim(2), points_per_dim);
    % create all combinations of grid points
    combs = combvec(xs, ys);
    
    % evaluate function at grid points
    zs = arrayfun(@(i) exercise_3_function(combs(:,i)), 1:size(combs, 2));
    
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
    s = @(x) sum(arrayfun(@(i) lambda(i).*phi_r(norm(x - combs(:,i))), 1:N));
    
    % create test points
    test_points = (lhsdesign(100, 2) * lim(2))';
    
    % evaluate surrogate model and true function at test points
    f_test = arrayfun(@(i) exercise_3_function(test_points(:,i)), 1:size(test_points, 2));
    s_test = arrayfun(@(i) s(test_points(:,i)), 1:size(test_points, 2));
    % calculate RMSE
    error = rmse(f_test, s_test);
    
    % plot
    figure;
    fmesh(@(x,y) s([x;y]), lim);
    xlim([-4 4])
    ylim([-4 4])
    hold on;
    scatter3(combs(1,:), combs(2,:), zs, 'k', 'filled');
    hold off
    title({['Surrogate model for N = ', num2str(N), ', RMSE = ', num2str(error)]});
    % saveas(gcf, ['d11e3_', num2str(N), '.png']);
    drawnow
end