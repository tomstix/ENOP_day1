clear; clc; close all;

load('exercise_1_data.mat');
x_values = linspace(0, 15, width(Y));

n = 8;

% create basis functions
basis_functions = cell(1, 2*n+1);
basis_functions{1} = @(x) 1;
for i = 2:2:2*n+1
    basis_functions{i} = @(x) sin(i/2*x);
    basis_functions{i+1} = @(x) cos(i/2*x);
end

% fill X matrix
X = zeros(length(x_values), 2*n+1);
for i = 1:length(x_values)
    for j = 1:2*n+1
        X(i, j) = basis_functions{j}(x_values(i));
    end
end

% solve linear system
lambda = (X'*X) \ (X'*Y');

% create surrogate model
s = @(x) sum(arrayfun(@(i) lambda(i).*basis_functions{i}(x), 1:2*n+1));

% plot
figure
plot(x_values, Y, 'bx');
hold on
fplot(@(x) s(x), [0 15], 'k');
title({['Regression model with n = ', num2str(n)]});
legend('True function', 'Surrogate model');
xlabel('x');
ylabel('f(x)');
