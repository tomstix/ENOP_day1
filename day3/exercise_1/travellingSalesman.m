% this variable is used to set the number of cities
numCities = 50;
% this variable is used to set the size of the grid
gridSize = 100;

% create random cities (coordinates) within the grid
cities = createRandomCities(numCities, gridSize, gridSize);

% generate the path and the length of the path
[path, length] = evaluateGraph(cities);

% create the initial plot
figure
p = plot(path(:, 1), path(:, 2), 'o-');
% store gca so we can change the title later
ha = gca;
xlim([0, gridSize])
ylim([0, gridSize])
title(['Path length: ', num2str(length)])
% set the data sources for the plot so we can update the plot later
p.XDataSource = 'path(:, 1)';
p.YDataSource = 'path(:, 2)';

% this is to break the loop if we are stuck
numUnsuccessfulSwaps = 0;
maxUnsuccessfulSwaps = 10000;

while numUnsuccessfulSwaps < maxUnsuccessfulSwaps
    % swap two random cities in the path
    newPath = randomSwap(path);
    [newPath, newLength] = evaluateGraph(newPath);
    
    % update the path if the new path is shorter
    if newLength < length
        path = newPath;
        length = newLength;
        % update the plot
        refreshdata
        drawnow
        title(['Path length: ', num2str(length)])
        numUnsuccessfulSwaps = 0;
        % pause(0.5)
    else
        numUnsuccessfulSwaps = numUnsuccessfulSwaps + 1;
    end
end

function newPath = randomSwap(oldPath)
% copy the old path to the new path
newPath = oldPath;

% generate two random indices that are not the first and the last one
i = randi([2, size(newPath, 1) - 1]);
j = randi([2, size(newPath, 1) - 1]);

% swap the cities at the two indices
newPath([i, j], :) = newPath([j, i], :);
end

function randomCities = createRandomCities(n, xMax, yMax)
randomCities = [randi(xMax, n, 1), randi(yMax, n, 1)];
end