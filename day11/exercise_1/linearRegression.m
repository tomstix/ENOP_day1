clear; clc; close all;

% Dataset
dataSet = importdata('exercise_1_data.mat');

% Initial plot of dataset
xAxis = linspace(0,15,size(dataSet,2));
% figure
% scatter(xAxis, dataSet, 50, 'x', 'LineWidth', 1.5, 'MarkerEdgeColor', [0 0 1]);
% hold on

% linear regression
n_max = 20;
results = zeros(n_max, size(xAxis,2));
errors = zeros(1, n_max);
for n = 1:n_max
    lambda = findLambda(n, xAxis, dataSet');
    results(n,:) = regressionModel(xAxis, n, lambda);
    errors(n) = rmse(dataSet, results(n,:));
    % plot(xAxis, results(n,:), 'LineWidth', 1);
end 

figure
subplot(2, 3, 1)
scatter(xAxis, dataSet, 50, 'x', 'LineWidth', 1.5, 'MarkerEdgeColor', [0 0 1]);
grid on
xlim([0 15])
ylim([-20 15])
title("Data Points")
subplot(2, 3, 2)
n = 2;
scatter(xAxis, dataSet, 50, 'x', 'LineWidth', 1.5, 'MarkerEdgeColor', [0 0 1]);
hold on
plot(xAxis, results(n,:), 'LineWidth', 1, 'Color', [0 0 0]);
grid on
xlim([0 15])
ylim([-20 15])
title({['Regression Model for {\it n} = ', num2str(n)]})
subplot(2, 3, 3)
n = 5;
scatter(xAxis, dataSet, 50, 'x', 'LineWidth', 1.5, 'MarkerEdgeColor', [0 0 1]);
hold on
plot(xAxis, results(n,:), 'LineWidth', 1, 'Color', [1 0 1]);
grid on
xlim([0 15])
ylim([-20 15])
title({['Regression Model for {\it n} = ', num2str(n)]})
subplot(2, 3, 4)
n = 8;
scatter(xAxis, dataSet, 50, 'x', 'LineWidth', 1.5, 'MarkerEdgeColor', [0 0 1]);
hold on
plot(xAxis, results(n,:), 'LineWidth', 1, 'Color', [1 0 0]);
grid on
xlim([0 15])
ylim([-20 15])
title({['Regression Model for {\it n} = ', num2str(n)]})
subplot(2, 3, 5)
n = 20;
scatter(xAxis, dataSet, 50, 'x', 'LineWidth', 1.5, 'MarkerEdgeColor', [0 0 1]);
hold on
plot(xAxis, results(n,:), 'LineWidth', 1, 'Color', [0 1 0]);
grid on
xlim([0 15])
ylim([-20 15])
title({['Regression Model for {\it n} = ', num2str(n)]})
subplot(2, 3, 6)
n = 20;
plot(1:1:size(errors,2), errors, 'LineWidth', 1, 'Color', [0 0 1]);
grid on
title("Modeling Error versus Model Order {\it n}")

function lambda = findLambda(n, x, y)
    X = zeros(size(y,1), n);
    for i = 1:size(X,1)
        X(i,1) = 1;
        for j = 2:2:(n*2)
            X(i,j) = sin((j/2)*x(i));
            X(i,j + 1) = cos((j/2)*x(i));
        end
    end
    lambda = ((X'*X)^-1*X'*y)';
end

function y = regressionModel(x, n, lambda)
    y = lambda(1) * 1;
    for i = 2:2:(n*2)
        y = y + lambda(i) .* sin((i/2).*x) + lambda(i + 1) .* cos((i/2).*x);
    end
end