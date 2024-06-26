clear; clc; close all;
% suppress fplot warning. didn't get this to be vectorized
warning('off', 'MATLAB:fplot:NotVectorized');

data = [
    0.0000,  1.0253;
    0.1000,  0.8702;
    0.2000,  0.5632;
    0.3000,  0.1260;
    0.4000, -0.2467;
    0.5000, -0.5407;
    0.6000, -0.6864;
    0.7000, -0.5969;
    0.8000, -0.4323;
    0.9000, -0.0802;
    1.0000,  0.2176;
    1.1000,  0.4952;
    1.2000,  0.6125;
    1.3000,  0.5570;
    1.4000,  0.4531;
    1.5000,  0.2293
    ];

cs = 2 + zeros(1, 100);
cs = cs .^ linspace(1, 10, 100);

% set up radial basis function
phi_r = @(r, c) exp(-c .* r.^2);

% group data
data_group_idx = (1:height(data)/2)';
data_group_idx = [data_group_idx; data_group_idx];
% data could also be grouped like this but it will create bigger gaps and might not be as good
% data_group_idx = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];

avg_errors = zeros(size(cs));
for c = cs
    % perform cross-validation
    errors = zeros(height(data)/2, 1);
    for k = 1:height(data)/2
        % use data from data_group k as test data
        test_data = data(data_group_idx == k, :);
        % the rest as training data
        train_data = data(data_group_idx ~= k, :);
        
        % do rbf interpolation
        lambda = rbf(train_data, @(x) phi_r(x, c));
        s = rbf_model(train_data, lambda, @(x)phi_r(x,c));
        
        % evaluate surrogate model and true function at test points
        f_test = test_data(:, 2);
        s_test = arrayfun(@(i) s(test_data(i, 1)), 1:size(test_data, 1))';
        
        % calculate RMSE
        errors(k) = rmse(f_test, s_test);
    end
    % calculate average error for cross-validation
    avg_errors(c == cs) = mean(errors);
    fprintf("c = %.2f, avg error = %.4f\n", c, avg_errors(c == cs));
end

% find best c
[~, best_c_idx] = min(avg_errors);
best_c = cs(best_c_idx);

% create model for best c
lambda = rbf(data, @(x) phi_r(x, best_c));
s_best = rbf_model(data, lambda, @(x)phi_r(x,best_c));
% create model for c = 4
lambda = rbf(data, @(x) phi_r(x, 4));
s_4 = rbf_model(data, lambda, @(x)phi_r(x,4));
% create model for c = 1024
lambda = rbf(data, @(x) phi_r(x, 1024));
s_1024 = rbf_model(data, lambda, @(x)phi_r(x,1024));

% plot average errors
figure;
plot(cs, avg_errors);
grid on
xscale log
yscale log
xlim("padded")
ylim("padded")
xlabel("c")
ylabel("Average RMSE")

% plot examples and best model
figure;
set(gcf, 'Position', [100, 100, 1200, 400]);

xli = [-0.1 1.6];
yli = [-1 1.5];

subplot(1, 3, 1);
hold on;
plot(data(:, 1), data(:, 2), 'o');
fplot(s_4, [0 1.5], 'r');
grid on;
hold off;
title("Model for c = 4");
xlim(xli)
ylim(yli)
xlabel("$x$", 'Interpreter', 'latex');
ylabel("$f(x)$", 'Interpreter', 'latex');

subplot(1, 3, 2);
hold on;
plot(data(:, 1), data(:, 2), 'o');
fplot(s_best, [0 1.5], 'r');
grid on;
hold off;
title(sprintf("Best model for c = %.2f", best_c));
xlim(xli)
ylim(yli)
xlabel("$x$", 'Interpreter', 'latex');
ylabel("$f(x)$", 'Interpreter', 'latex');

subplot(1, 3, 3);
hold on;
plot(data(:, 1), data(:, 2), 'o');
fplot(s_1024, [0 1.5], 'r');
grid on;
hold off;
title("Model for c = 1024");
xlim(xli)
ylim(yli)
xlabel("$x$", 'Interpreter', 'latex');
ylabel("$f(x)$", 'Interpreter', 'latex');

function lambda = rbf(data, phi_r)
% function to calculate lambda for radial basis function
% create Phi matrix
N = size(data, 1);
Phi = zeros(N, N);
for i = 1:N
    for j = 1:N
        Phi(i, j) = phi_r(norm(data(i, 1) - data(j, 1)));
    end
end

% calculate lambda
lambda = Phi \ data(:, 2);
end

function y = rbf_model(data, lambda, phi_r)
% function to create model for rbf from data, lambda and phi_r
y = @(x) sum(arrayfun(@(i) lambda(i).*phi_r(norm(x - data(i, 1))), 1:size(data, 1)));
end
