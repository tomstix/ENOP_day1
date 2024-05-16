clear; clc; close all;

% Control variables
x_1 = 1;
xAxis = -5:0.1:5;

% High- and low-fidelity function definitions
func = @(x) exp(x ./ 3) .* 0.1 .* x.^2 ./ (1 + 0.1 .* x.^2);
cunc = @(x) exp(x ./ 3);

% Initial plot of high and low fidelity model
% figure
% plot(xAxis, func(xAxis), "Color", [0 0 1], "LineWidth", 1);
% hold on
% plot(xAxis, cunc(xAxis), "Color", [0 0 0], "LineWidth", 1);
% legend('High-fidelity model', 'Low-fidelity model');
% grid on
% xticks(xAxis(1):1:xAxis(end));
% xlim([xAxis(1) xAxis(end)]);
% ylim([-2 6]);

% Tangent function and gradient at x^(1)
[func_, func_grad_x_1] = firstOrderTaylor(func, x_1, 0.1);
[cunc_, cunc_grad_x_1] = firstOrderTaylor(cunc, x_1, 0.1);

% Beta-correlation method
aunc = @(x) ((func(x_1) / cunc(x_1)) + ((func_grad_x_1 * cunc(x_1) - cunc_grad_x_1 * func(x_1))/cunc(x_1)^2) * (x - x_1));
sunc = @(x) aunc(x).*cunc(x);

% Plot
figure
subplot(1, 2, 1)
plot(xAxis, func(xAxis), "Color", [0 0 1], "LineWidth", 1);
hold on
plot(xAxis, cunc(xAxis), "Color", [0 0 0], "LineWidth", 1);
scatter(x_1, func(x_1), 25, 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1)
plot(xAxis, sunc(xAxis), "Color", [1 0 0], "LineWidth", 1);
plot(xAxis, func_(xAxis), "Color", [1 0 0], "LineWidth", 1, "LineStyle", "--");
legend('High-fidelity model {\it f(x)}', 'Low-fidelity model {\it c(x)}', '{\it f(x^{(1)})}', 'Surrogate model {\it s(x)}', 'Tangent {\it f(x^{(1)}) + f''(x-x^{(1)})}', 'Location','northwest');
grid on
xticks(xAxis(1):1:xAxis(end));
xlim([xAxis(1) xAxis(end)]);
ylim([-2 6]);
title(['Beta-correlation method, {\it x^{(1)} = ', num2str(x_1),'}']);

% General response correction
sunc = @(x) func(x_1) + (func_grad_x_1 / cunc_grad_x_1)*(cunc(x) - cunc(x_1));

% Plot
subplot(1, 2, 2)
plot(xAxis, func(xAxis), "Color", [0 0 1], "LineWidth", 1);
hold on
plot(xAxis, cunc(xAxis), "Color", [0 0 0], "LineWidth", 1);
scatter(x_1, func(x_1), 25, 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1)
plot(xAxis, sunc(xAxis), "Color", [1 0 0], "LineWidth", 1);
plot(xAxis, func_(xAxis), "Color", [1 0 0], "LineWidth", 1, "LineStyle", "--");
legend('High-fidelity model {\it f(x)}', 'Low-fidelity model {\it c(x)}', '{\it f(x^{(1)})}', 'Surrogate model {\it s(x)}', 'Tangent {\it f(x^{(1)}) + f''(x-x^{(1)})}', 'Location','northwest');
grid on
xticks(xAxis(1):1:xAxis(end));
xlim([xAxis(1) xAxis(end)]);
ylim([-2 6]);
title(['General response correction, {\it x^{(1)} = ', num2str(x_1),'}']);
