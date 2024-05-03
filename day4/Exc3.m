
result1 = trust_region_method(@quadratic_function, (1:10)')
disp(result1);

result2 = trust_region_method(@sin, 2);
disp(result2);

result3 = trust_region_method(@expon_function, [5,3]')
disp(result3);

result4 = trust_region_method(@rosenbr_function, [0,0]')
disp(result4)

% Prepare subplot
figure;

% % Plot result for quadratic function
% xquad_range=linspace(-5,5,100);
% subplot(2, 2, 2);
% n = length(result2);
% plot(xquad_range,quadratic_function(xquad_range),'r-');
% hold on;
% plot(result2, quadratic_function(result2), 'bo');
% title('quadratic');
% xlabel('x');
% ylabel('f(x)');
% grid on;
% hold off;

% Plot result for quadratic function
subplot(2, 2, 1);
xquad_range = linspace(-10, 10, 100);
yquad_values = arrayfun(@quadratic_function, xquad_range);
plot(xquad_range, yquad_values, 'r-'); % Function curve
hold on;
plot((1:10)', quadratic_function((1:10)'), 'bo'); % Start point
plot(result1, quadratic_function(result1), 'go'); % End point
title('quadratic');
xlabel('x');
ylabel('f(x)');
legend('Function curve', 'Start point', 'End point');
grid on;
hold off;
% % Plot result for sin(x)
% xsin_range=linspace(3.5,6,10);
% subplot(2, 2, 1);
% plot(xsin_range,sin(xsin_range),'r-');
% hold on;
% plot(result1, sin(result1), 'bo');
% title('sin(x)');
% xlabel('x');
% ylabel('sin(x)');
% grid on;
% hold off;

% Plot result for sin(x)
xsin_range = linspace(3.5, 6, 100);
subplot(2, 2, 2);
plot(xsin_range, sin(xsin_range), 'r-'); % Function curve
hold on;
plot(2, sin(2), 'bo'); % Start point
plot(result2, sin(result2), 'go'); % End point
title('sin(x)');
xlabel('x');
ylabel('sin(x)');
legend('Function curve', 'Start point', 'End point');
grid on;
hold off;

% % Plot result for exponential function
% xexp_range=linspace(-2,2,10);
% yexp_range=arrayfun(@(x) expon_function([x, 3]), xexp_range);
% subplot(2, 2, 3);
% plot(xexp_range, yexp_range, 'r-');
% hold on;
% plot(result3(1), result3(2), 'bo');
% title('exp function');
% xlabel('x');
% ylabel('h(x)');
% grid on;

% Plot result for exponential function
subplot(2, 2, 3);
xexp_range = linspace(-10, 10, 100);
yexp_values = arrayfun(@(x) expon_function([x, 3]), xexp_range);
plot(xexp_range, yexp_values, 'r-'); % Function curve
hold on;
plot(5, expon_function([5, 3]), 'bo'); % Start point
plot(result3(1), expon_function(result3), 'go'); % End point
title('exp function');
xlabel('x');
ylabel('h(x)');
legend('Function curve', 'Start point', 'End point');
grid on;
hold off;

% % Plot result for Rosenbrock function
% subplot(2, 2, 4);
% plot(result4(1), result4(2), 'bo');
% title('Rosenbrock');
% xlabel('x');
% ylabel('r(x)');
% grid on;

% Plot result for Rosenbrock function
subplot(2, 2, 4);
xrosen_range = linspace(-2, 2, 100);
yrosen_values = arrayfun(@(x) rosenbr_function([x, x]), xrosen_range);
plot(xrosen_range, yrosen_values, 'r-'); % Function curve
hold on;
plot(0, rosenbr_function([0, 0]), 'bo'); % Start point
plot(result4(1), rosenbr_function(result4), 'go'); % End point
title('Rosenbrock');
xlabel('x');
ylabel('r(x)');
legend('Function curve', 'Start point', 'End point');
grid on;
hold off;


function y = quadratic_function(x);
    y= sum(x.^2);
end

function h= expon_function(x);
    h= exp(x(1)/5+x(2)/2)+x(1)^2+x(2)^2;
end

function r= rosenbr_function(x);
    r=(1-x(1)).^2+100*(x(2)-x(1).^2).^2;
end

function x = trust_region_method(g, x0, delta_max, delta_min, tol, max_iter, eta)
    if nargin < 2, x0 = 2; end
    if nargin < 3, delta_max = 1.0; end
    if nargin < 4, delta_min = 1e-6; end
    if nargin < 5, tol = 1e-6; end
    if nargin < 6, max_iter = 100; end
    if nargin < 7, eta = 0.2; end

    x = x0;
    delta = delta_max;
    for k = 1:max_iter
        grad = finite_difference_gradient(g, x);
        if abs(grad) < tol
            break;
        end

        linear_model = @(s) g(x) + grad * s;
        s_star = golden_section_search(linear_model, -delta, delta);
        x_new = x + s_star;
        actual_reduction = g(x) - g(x_new);
        predicted_reduction = linear_model(0) - linear_model(s_star);
        if predicted_reduction ~= 0
            rho = actual_reduction / predicted_reduction;
        else
            rho = 0;
        end

        if rho < 0.25
            delta = max(0.25 * delta, delta_min);
        elseif rho > 0.75 && abs(s_star) == delta
            delta = min(2.0 * delta, delta_max);
        end

        if rho > eta
            x = x_new;
        end
    end
end

function grad = finite_difference_gradient(g, x)
    epsilon = 1e-8;
    grad = (g(x + epsilon) - g(x)) / epsilon;
end

function x = golden_section_search(f, a, b, tol)
    if nargin < 4, tol = 1e-5; end
    gr = (sqrt(5) + 1) / 2;
    c = b - (b - a) / gr;
    d = a + (b - a) / gr;
    while abs(c - d) > tol
        if f(c) < f(d)
            b = d;
        else
            a = c;
        end
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
    end
    x = (b + a) / 2;
end