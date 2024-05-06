function [newGridSize, bestNeighbour] = searchBestNeighbour(points, fn, gridSize, ax)
current = points(end,:);
newGridSize = gridSize;
bestNeighbour = current;
foundBetter = false;
while ~foundBetter
    % create all neighbours
    neighbours(1,:) = [current(1) + newGridSize, current(2)];
    neighbours(2,:) = [current(1) - newGridSize, current(2)];
    neighbours(3,:) = [current(1), current(2) + newGridSize];
    neighbours(4,:) = [current(1), current(2) - newGridSize];
    % check all neighbours
    for i = 1:4
        % plot the line between the current point and the neighbour
        pl = plot(ax, [current(1), neighbours(i,1)], [current(2), neighbours(i,2)], 'g', 'LineWidth', 1);
        if fn(neighbours(i,1), neighbours(i,2)) < fn(bestNeighbour(1), bestNeighbour(2))
            bestNeighbour = neighbours(i,:);
            foundBetter = true;
        end
        % pause(0.05);
        delete(pl); % delete the line again
    end
    if ~foundBetter
        % if no better neighbour was found, reduce the grid size
        newGridSize = newGridSize / 3;
    end
end