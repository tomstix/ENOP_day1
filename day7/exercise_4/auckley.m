function f = auckley(x)
    if height(x) > 3
        error('This function is not implemented for more than 3D input')
    end
    switch height(x)
        case 1
            f = -20*exp(-0.2*sqrt(x(1,:).^2)) - exp(cos(2*pi*x(1,:))) + 21;
        case 2
            f = -20*exp(-0.2*sqrt(1/2 * (x(1,:).^2 + x(2,:).^2))) - exp(1/2 * (cos(2*pi*x(1,:)) + cos(2*pi*x(2,:)))) + 21;
        case 3
            f = -20*exp(-0.2*sqrt(1/3 * (x(1,:).^2 + x(2,:).^2 + x(3,:).^2))) - exp(1/3 * (cos(2*pi*x(1,:)) + cos(2*pi*x(2,:)) + cos(2*pi*x(3,:)))) + 21;
    end
end