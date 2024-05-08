clear;

data = [...
    1,2.3743;
    2,1.1497;
    3,0.7317;
    4,0.5556;
    5,0.4675;
    6,0.4157;
    7,0.3807;
    8,0.3546;
    9,0.3337;
    10,0.3164
    ];

t = data(:,1);
y = data(:,2);
fun = @(x) x(1) .* exp(-x(2) * t) + x(3) .* exp(-x(4) * t) - y;

lb = zeros(1,4);
ub = ones(1,4) * 10;

x0 = ones(1,4);
y0 = fun(x0) + y;

x = lsqnonlin(fun, x0, lb, ub);

figure
hold on
grid on
plot(t, y, 'k*');
plot(t, y0, 'bo');
plot(t, fun(x)+y, 'ro');