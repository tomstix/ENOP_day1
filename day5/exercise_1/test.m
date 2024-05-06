clear;

eps = 1e-3;

selected_function = 2;
switch selected_function
    case 1
        fn = @(x,y) 2*x.^2 + 3*y.^2 - 3.*x.*y + x;
        x = [5, 8];
        lim = [-2, 8];
        exactMin = [0, 0];
    case 2
        fn = @(x,y) (1-x).^2 + 5*(x-y.^2).^2;
        x = [0, 0];
        lim = [-0.5, 1.5];
        exactMin = [1, 1];
    case 3
        fn = @(x,y) (x+2.*y) .* exp(1-0.9.*exp(-0.3.*(x-2.5).^2 - 2.*(y-3.5).^2)) .* (1-0.9.*exp(-(x-3).^2 - (y-3).^2));
        x = [4, 2];
        lim = [1, 5];
        exactMin = [3, 3];
    case 4
        fn = @(x,y) exp(x/5) + exp(y/3);
        x = [5, 8];
        lim = [-10 10];
end

figure
contour = fcontour(fn, lim, "LevelStep", 2);
xlim(lim)
ylim(lim)
hold on
ax = gca;
drawnow

points = x;
visualizePath(ax, points, false);
ax = gca;
gridSize = 1;
iteration_count = 0;
while gridSize > eps
    iteration_count = iteration_count + 1;
    [gridSize, bestNeighbour] = searchBestNeighbour(points, fn, gridSize, ax);
    x = bestNeighbour;
    points = [points; x];
    % check limits
    if x(1) < lim(1) || x(1) > lim(2) || x(2) < lim(1) || x(2) > lim(2)
        break;
    end
    visualizePath(ax, points, false);
    fprintf("Iteration %d: x = %d, y = %d, f = %d, gridSize = %d\n", iteration_count, x(1), x(2), fn(x(1), x(2)), gridSize);
    pause(0.1);
end
visualizePath(ax, points, true);
