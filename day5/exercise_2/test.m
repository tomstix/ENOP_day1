clear;

eps = 1e-3;

selected_function = 1;
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

ax = contourPlot(fn, lim);

points = x;
visualizePath(ax, points, false);

gridSize = 1;
iteration_count = 0;
while gridSize > eps
    iteration_count = iteration_count + 1;
    grad = getPseudoGradient(fn, points(end,:), gridSize);
    q = quiver(ax, points(end,1), points(end,2), -grad(1), -grad(2), gridSize * 0.8, 'color', 'r', 'LineWidth', 1);
    xNew = lineSearch(fn, points(end,:), grad, gridSize);
    if fn(xNew(1), xNew(2)) < fn(points(end,1), points(end,2))
        bestNeighbour = xNew;
    else % fall back to neighbour search
        [gridSize, bestNeighbour] = searchBestNeighbour(points, fn, gridSize, ax);
    end
    x = bestNeighbour;
    points = [points; x];
    % check limits
    if x(1) < lim(1) || x(1) > lim(2) || x(2) < lim(1) || x(2) > lim(2)
        break;
    end
    visualizePath(ax, points, false);
    fprintf("Iteration %d: x = %d, y = %d, f = %d, gridSize = %d\n", iteration_count, x(1), x(2), fn(x(1), x(2)), gridSize);
    pause(0.1);
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

function point = lineSearch(fn, point, grad, gridSize)
% store the original point to correctly snap to grid
orig = point;
while true
    nextPoint = point - grad * gridSize;
    % go on until we find a point where the function value increases
    if fn(nextPoint(1), nextPoint(2)) > fn(point(1), point(2))
        break;
    end
    % this point is better, so we store it
    point = nextPoint;
end
% snap the final point to the grid
point = snapToGrid(orig, point, gridSize);
end