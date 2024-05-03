% Variable that decides which function Newton's method is used on
functionSelect = 4;
% Change the above variable and rerun the program to see Newton's method
% used on the 4 different functions.
% Valid values are 1, 2, 3 and 4, which function is which can be seen below

% Defining the functions
funck{1} = @(x) x.^2;
diffFunck{1} = @(x) 2*x;
x_0{1} = 5;
funck{2} = @(x) log(x);
diffFunck{2} = @(x) 1/x;
x_0{2} = 0.1;
funck{3} = @(x) x.^4;
diffFunck{3} = @(x) 4*x^3;
x_0{3} = 5;
funck{4} = @(x) x.^(0.5)-2;
diffFunck{4} = @(x) 1/(2*sqrt(x));
x_0{4} = 10;

% Stepsize for Finite difference derivatives
h = 0.0001;

% Setting up plot
x = x_0{functionSelect};
x_axis = -1:0.1:x_0{functionSelect}+2.5;
plot(x_axis,funck{functionSelect}(x_axis), 'LineWidth', 2)
hold on
plot(x_axis, 0*x_axis, 'LineWidth', 1, 'Color', [0 1 0])
grid on
xlim([x_axis(1) x_axis(end)])
ylim([funck{functionSelect}(0)-(funck{functionSelect}(x_axis(end))*0.1) funck{functionSelect}(x_axis(end))+(funck{functionSelect}(x_axis(end))*0.1)])
title(sprintf('Iteration: %d   x=%d    f=%d', 0, x_n, funck{functionSelect}(x_n)));
drawnow

% Newtons sequence
Iterations = 0;
while funck{functionSelect}(x) ~= 0
    x_n = x-funck{functionSelect}(x)/diffFunck{functionSelect}(x);

    % Plotting the sequence
    plot([x_n x], [0 funck{functionSelect}(x)], 'LineWidth', 1, 'Color', [1 0 0])
    plot([x_n x_n], [0 funck{functionSelect}(x_n)], 'LineWidth', 0.5, 'Color', [1 0 0])
    scatter([x_n x x_n], [0 funck{functionSelect}(x) funck{functionSelect}(x_n)], 5, [1 0 0], 'filled')
    Iterations = Iterations + 1;
    title(sprintf('Iteration: %d   x=%d    f=%d', Iterations, x_n, funck{functionSelect}(x_n)));
    drawnow limitrate

    x = x_n;
end
% Plot final x
scatter(x, 0, 36, [1 0 0], 'filled')