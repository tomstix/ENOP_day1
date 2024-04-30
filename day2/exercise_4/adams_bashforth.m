function [x, y] = adams_bashforth(f, x0, y0, h, x)
    n = length(x);
    y = zeros(1, n);
    x(1) = x0;
    y(1) = y0;
    % modified euler method
    y(2) = y(1) + h * f(x(1), y(1));
    y(3) = y(2) + h * f(x(2), y(2));
    for i = 2:n-2
        y(i+2) = y(i+1) + 1.5 * h * f(x(i+1), y(i+1)) - 0.5 * h * f(x(i), y(i));
    end
end