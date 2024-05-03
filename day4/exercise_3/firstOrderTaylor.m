function model = firstOrderTaylor(f, x0, h)
if width(x0) > 1
    x0 = x0';
    if width(x0) > 1
        error('x0 must be a vector');
    end
end
dim = length(x0);
grad = zeros(dim, 1);
for i = 1:dim
    m = zeros(dim, 1);
    m(i) = 1;
    fwd = x0 + m*h;
    bwd = x0 - m*h;
    grad(i) = (f(fwd) - f(bwd))/(2*h);
end
model = @(x) f(x0) + grad'*(x-x0);
end