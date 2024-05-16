function [model, grad] = firstOrderTaylor(f, x0, h)
% Check if x0 is a vector
if width(x0) > 1
    x0 = x0';
    if width(x0) > 1
        error('x0 must be a vector');
    end
end
% Calculate the gradient using central differences
dim = length(x0);
grad = zeros(dim, 1);
for i = 1:dim
    m = zeros(dim, 1);
    m(i) = 1;
    fwd = x0 + m*h;
    bwd = x0 - m*h;
    grad(i) = (f(fwd) - f(bwd))/(2*h);
end
% Return an anonymous function that represents the first order Taylor model
model = @(x) f(x0) + grad'*(x-x0);
end