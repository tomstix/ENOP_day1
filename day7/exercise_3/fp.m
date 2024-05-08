function f = fp(x)
    if height(x) > 2
        error('This function is not implemented for more than 2D input')
    end
    norms = sqrt(sum(x.^2));
    f = -exp(-(norms.^2)./2);
    for i = 1:height(x)
        f = f .* cos(10*pi.*x(i,:));
    end
end