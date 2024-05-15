clear; clc; close all;

% Dataset
dataSet(:,1) = (0:0.2:10)';
dataSet(:,2) = [2.9004; 2.8718; 2.6121; 2.5500; 2.3605; 2.0048; 1.8463; 1.5930; 1.2475; 1.1892; 1.0805; 0.9076; 0.7522; 0.7779; 0.6789; 0.6358; 0.5275; 0.5860; 0.6809; 0.7591; 0.7995; 0.8182; 0.9791; 0.9631; 1.0600; 1.1088; 1.1188; 1.0386; 0.9028; 0.8256; 0.6602; 0.6062; 0.4935; 0.3788; 0.2423; 0.1860; 0.1158; 0.1396; 0.1260; 0.1131; 0.0669; 0.0647; 0.0174; 0.0864; 0.0424; 0.0190; 0.0805; 0.0670; 0.0761; 0.0242; 0.0561];

% Initial plot of dataset
figure
scatter(dataSet(:,1),dataSet(:,2), 'x', 'LineWidth', 1.5, 'MarkerEdgeColor', [0 0 0]);
hold on

% Function definition
parameters = [1, 2, 1.5, -5];
% funck = @(a, b, c, d, x) (cos(a*x) + b) .* (x+c).^d;
funck = @(a, b, c, d, x) (cos(a*x) + b) .* (1./(1+exp(c.*(x+d))));

% Find function that takes a, b, c and d as parameters and returns the
% difference between the model and the actual data
d = dataSet(:,1);
y = dataSet(:,2);
fun = @(x) funck(x(1), x(2), x(3), x(4), d) - y;

% Optimize parameters
parameters = lsqnonlin(fun, parameters);

% Plot setup
plot(d, funck(1, 2, 1.5, -5, d), 'LineWidth', 1.5, 'Color', [0 0 1]);
plot(d, funck(parameters(1), parameters(2), parameters(3), parameters(4), d), 'LineWidth', 1.5, 'Color', [1 0 0]);
legend('Data', 'Approximation function (initial)', 'Approximation function (optimized)', 'Location','northeast');
grid on
xlabel('{\it x}') 
ylabel('{\it y}')