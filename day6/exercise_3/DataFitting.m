clear
clc

% Dataset
data = [1 2.3743; 2 1.1497; 3 0.7317; 4 0.5556; 5 0.4675; 6 0.4157; 7 0.3807; 8 0.3546; 9 0.3337; 10 0.3164];
bounds = [0 10];

% Function definition
parameters = [1 1 1 1];
funck = @(a, b, c, d, t) a .* exp(-b * t) + c .* exp(-d * t);

% Calculate initial approximation values
d = data(:,1);
y = data(:,2);

% Find function that takes a, b, c and d as parameters and returns the
% difference between the model and the actual data
fun = @(x) funck(x(1), x(2), x(3), x(4), d) - y;

% Optimize parameters
lb = ones(1, size(parameters, 2)) * bounds(1);
ub = ones(1, size(parameters, 2)) * bounds(2);
parameters = lsqnonlin(fun, parameters, lb, ub);

% Plot setup
scatter(d, data(:,2), 'x', 'LineWidth', 1.5, 'MarkerEdgeColor', [0 0 0]);
hold on
scatter(d, funck(1, 1, 1, 1, d), 'LineWidth', 1.5, 'MarkerEdgeColor', [0 0 1]);
scatter(d, funck(parameters(1), parameters(2), parameters(3), parameters(4), d), 'LineWidth', 1.5, 'MarkerEdgeColor', [1 0 0]);
legend('Data', 'Approximation function (initial)', 'Approximation function (optimized)', 'Location','northeast');
grid on
xlabel('{\it t}') 
ylabel('{\it y}')