clear;

% select the break condition
eps = 1e-3;
max_iterations = 100;

% select the initial grid size
gridSize = 1;

% select which function to use
selected_function = 1;
switch selected_function
    case 1
        fn = @(x,y) 2*x.^2 + 3*y.^2 - 3.*x.*y + x;
        x = [5, 8]; % starting point
        lim = [-2, 8]; % limits
        viewSettings = [-37.5, 50]; % view settings for 3D plot
    case 2
        fn = @(x,y) (1-x).^2 + 5*(x-y.^2).^2;
        x = [0, 0];
        lim = [-0.5, 1.5];
        viewSettings = [-60, 45];
    case 3
        fn = @(x,y) (x+2.*y) .* exp(1-0.9.*exp(-0.3.*(x-2.5).^2 - 2.*(y-3.5).^2)) .* (1-0.9.*exp(-(x-3).^2 - (y-3).^2));
        x = [4, 2];
        lim = [1, 5];
        viewSettings = [-45, 65];
    case 4
        fn = @(x,y) exp(x/5) + exp(y/3);
        x = [5, 8];
        viewSettings = [-37.5, 30];
        lim = [-10, 10];
end

% create initial contour plot
figure
contour = fcontour(fn, lim, "LevelStep", 2);
xlim(lim)
ylim(lim)
hold on
ax = gca;
drawnow

points = x;
visualizePath(ax, points, false);
iteration_count = 0;
while gridSize > eps && iteration_count < max_iterations
    iteration_count = iteration_count + 1;
    fprintf("Iteration %d: x = %d, y = %d, f = %d, gridSize = %d\n", iteration_count, x(1), x(2), fn(x(1), x(2)), gridSize);
    [gridSize, bestNeighbour] = searchBestNeighbour(points, fn, gridSize, ax);
    % check limits
    if bestNeighbour(1) < lim(1) || bestNeighbour(1) > lim(2) || bestNeighbour(2) < lim(1) || bestNeighbour(2) > lim(2)
        break;
    end
    x = bestNeighbour;
    points = [points; x];
    visualizePath(ax, points, false);
    % pause(0.1);
end
visualizePath(ax, points, true);
