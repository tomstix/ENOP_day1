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

figure
subplot(1,3,1)
contour = fcontour(fn, lim, "LevelStep", 2);
xlim(lim)
ylim(lim)
hold on
ax2D = gca;
drawnow

subplot(1,3,2)
ax3D = gca;
fsurf(fn, lim);
hold on

subplot(1,3,3)
axZoom = gca;
fmesh(fn, lim);
hold on

points = x;
visualizePath(ax2D, points, false);
visualizePath3D(ax3D, points, fn, false);
visualizePath3D(axZoom, points, fn, false);
gridSize = 1;
iteration_count = 0;
while gridSize > eps
    iteration_count = iteration_count + 1;
    [gridSize, bestNeighbour] = searchBestNeighbour(points, fn, gridSize, ax2D);
    x = bestNeighbour;
    points = [points; x];
    % check limits
    if x(1) < lim(1) || x(1) > lim(2) || x(2) < lim(1) || x(2) > lim(2)
        break;
    end
    visualizePath(ax2D, points, false);
    visualizePath3D(ax3D, points, fn, false);
    visualizePath3D(axZoom, points, fn, false);
    
    if height(points) > 5
        lookAtPoints = points(end-5:end, :);
        lookAtPoints(:,3) = fn(lookAtPoints(:,1), lookAtPoints(:,2));
        xRange = min(1, max(lookAtPoints(:,1)) - min(lookAtPoints(:,1)));
        yRange = min(1, max(lookAtPoints(:,2)) - min(lookAtPoints(:,2)));
        zRange = min(1, max(lookAtPoints(:,3)) - min(lookAtPoints(:,3)));
        xlim(axZoom, [min(lookAtPoints(:,1)) - xRange, max(lookAtPoints(:,1)) + xRange]);
        ylim(axZoom, [min(lookAtPoints(:,2)) - yRange, max(lookAtPoints(:,2)) + yRange]);
        zlim(axZoom, [min(lookAtPoints(:,3)) - zRange, max(lookAtPoints(:,3)) + zRange]);
    end

    fprintf("Iteration %d: x = %d, y = %d, f = %d, gridSize = %d\n", iteration_count, x(1), x(2), fn(x(1), x(2)), gridSize);
    pause(0.1);
end
visualizePath(ax2D, points, true);
visualizePath3D(ax3D, points, fn, true);
visualizePath3D(axZoom, points, fn, true);
