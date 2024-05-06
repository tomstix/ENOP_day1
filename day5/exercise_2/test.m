clear;

% select the break condition
eps = 1e-3;
max_iterations = 100;

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

ax = contourPlot(fn, lim);

points = x;
visualizePath(ax, points, false);

gridSize = 1;
iteration_count = 0;
while gridSize > eps && iteration_count < max_iterations
    iteration_count = iteration_count + 1;
    grad = getPseudoGradient(fn, points(end,:), gridSize);
    q = quiver(ax, points(end,1), points(end,2), -grad(1), -grad(2), gridSize, 'color', 'r', 'LineWidth', 1);
    xNew = lineSearch(fn, points(end,:), grad, gridSize, lim);
    if fn(xNew(1), xNew(2)) < fn(points(end,1), points(end,2))
        bestNeighbour = xNew;
    else % fall back to neighbour search
        [gridSize, bestNeighbour] = searchBestNeighbour(points, fn, gridSize, ax);
    end
    % check limits
    if bestNeighbour(1) < lim(1) || bestNeighbour(1) > lim(2) || bestNeighbour(2) < lim(1) || bestNeighbour(2) > lim(2)
        break;
    end
    x = bestNeighbour;
    points = [points; x];
    visualizePath(ax, points, false);
    fprintf("Iteration %d: x = %d, y = %d, f = %d, gridSize = %d\n", iteration_count, x(1), x(2), fn(x(1), x(2)), gridSize);
    % pause(0.1);
    delete(q);
end
visualizePath(ax, points, true);

function ax = contourPlot(fn, lim)
figure
fcontour(fn, lim, "LevelStep", 2);
xlim(lim)
ylim(lim)
hold on
ax = gca;
drawnow
end

function newPoint = snapToGrid(pointOnGrid, pointToSnap, gridSize)
% snap a point to the grid given by pointOnGrid and gridSize
newPoint = pointOnGrid + round((pointToSnap - pointOnGrid) / gridSize) * gridSize;
end

function grad = getPseudoGradient(fn, point, gridSize)
% calculate a pseudo gradient in the grid using central differences
grad = zeros(1,2);
grad(1) = (fn(point(1) + gridSize, point(2)) - fn(point(1) - gridSize, point(2))) / (2 * gridSize);
grad(2) = (fn(point(1), point(2) + gridSize) - fn(point(1), point(2) - gridSize)) / (2 * gridSize);
% normalize
grad = grad / norm(grad);
end

function point = lineSearch(fn, point, grad, gridSize, lim)
% store the original point to correctly snap to grid
orig = point;
while true
    nextPoint = point - grad * gridSize;
    % go on until we find a point where the function value increases or we reach the limits
    if fn(nextPoint(1), nextPoint(2)) > fn(point(1), point(2)) || ...
            nextPoint(1) < lim(1) || nextPoint(1) > lim(2) || ...
            nextPoint(2) < lim(1) || nextPoint(2) > lim(2)
        break;
    end
    % this point is better, so we store it
    point = nextPoint;
end
% snap the final point to the grid
point = snapToGrid(orig, point, gridSize);
end