function f = rosenbrock(x)
if height(x)~= 2
    error('Rosenbrock function is only defined for 2D input (column vector)')
end
f = (1 - x(1,:).^2) + 100.*(x(2,:) - x(1,:).^2).^2;
end