function xmin = linearSearch(f, x, delta)
dim = length(x);
xmin = zeros(dim, 1);
for i = 1:dim
    x1 = x;
    x2 = x;
    x1(i) = x(i) - delta;
    x2(i) = x(i) + delta;
    f1 = f(x1);
    f2 = f(x2);
    if f1 < f2
        xmin(i) = x1(i);
    else
        xmin(i) = x2(i);
    end
end