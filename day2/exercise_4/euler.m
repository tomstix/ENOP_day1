function [x, y] = euler(f, x0, y0, h, x)
    n = length(x);
    y = zeros(1, n);
    x(1) = x0;
    y(1) = y0;
    for i = 1:n-1
        y(i+1) = y(i) + h * f(x(i), y(i));
    end
end