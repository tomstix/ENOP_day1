function f = ff(x)
    f(1, :) = 1 - exp( - sum( (x - 1/sqrt(8)).^2, 1) );
    f(2, :) = 1 - exp( - sum( (x + 1/sqrt(8)).^2, 1) );
end