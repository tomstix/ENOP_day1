%% Uncomment the function you want to test
% f = @(x, y) -y + 3*cos(3*x)*exp(-x);
% sol = @(x)sin(3*x).*exp(-x);
% x0 = 0;
% y0 = 0;

%%
f = @(x, y) y;
sol = @(x)exp(x);
x0 = 0;
y0 = 1;

%%
% initial conditions
% set step size here
hs = [0.2, 0.1, 0.05, 0.01];
for i = 1:length(hs)
    h = hs(i);
    % set range of x values
    x_values = 0:h:5;
    
    % solve using euler method
    [x, y] = euler(f, x0, y0, h, x_values);
    % solve using adams bashforth method
    [x2, y2] = adams_bashforth(f, x0, y0, h, x_values);
    
    subplot(2,2,i)
    hold on
    plot(x, y, "Color","black", "LineStyle", "--")
    plot(x2, y2, "Color","black", "LineStyle", "-.")
    fplot(sol, [min(x_values), max(x_values)], "Color","blue")
    lgd = legend("Euler", "Adams-Bashforth", "Analytical");
    lgd.Title.String = sprintf("h = %.2f", h);
    % saveas(gcf, sprintf("output/exercise_4_2_%d.png", i))
end