clear; clc;

n = 10; % select n here

x0 = ones(n,1) * 1/n; % starting point

% equality constraint
Aeq = ones(1,n);
beq = 1;

% lower bound
lb = zeros(n,1);

options = optimset('Display','iter'); % display output
x = fminimax(@black_box, x0, [], [], Aeq, beq, lb, [], [], options);

figure
set(gcf, 'Position', [100, 100, 400, 800])
subplot(2,1,1)
hold on;
% offset the x values for better visualization
xs = 1:n;
xs1 = xs - 0.1;
xs2 = xs + 0.1;
% plot the function arguments
plot(xs1, x0, 'ko');
plot(xs2, x, 'ro');
for i = 1:n
    line([i+0.1 i+0.1], [0 x(i)], 'Color', 'red');
    line([i-0.1 i-0.1], [0 x0(i)], 'Color', 'black');
end
grid on
xlim([1-0.2 n+0.2])
title('Function arguments')

subplot(2,1,2)
hold on;
% offset the x values for better visualization
xs = 1:n;
xs1 = xs - 0.1;
xs2 = xs + 0.1;
% plot the function arguments
values0 = black_box(x0);
values = black_box(x);
plot(xs1, values0, 'ko');
plot(xs2, values, 'ro');
for i = 1:n
    line([i+0.1 i+0.1], [0 values(i)], 'Color', 'red');
    line([i-0.1 i-0.1], [0 values0(i)], 'Color', 'black');
end
grid on
xlim([1-0.2 n+0.2])
title('Function values')