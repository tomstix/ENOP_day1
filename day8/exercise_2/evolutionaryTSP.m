clear
clc
close all

% this variable is used to set the number of cities
numCities = 50;
% this variable is used to set the size of the grid
gridSize = 100;

% create random cities (coordinates) within the grid
cities = [randi(gridSize, numCities, 1), randi(gridSize, numCities, 1)];

% generate the path and the length of the path
[path, length] = evaluateGraph(cities);
initialPath = path;
initialLength = length;

% Algorithm setup
maxIterations = 1000;
iterations = 1;
N = 100; % Population size
P = cell(N, 2); % Population
for i = 1:N
    cityOrder = randperm(numCities);
    cityOrder(find(cityOrder == 1, 1)) = [];
    individualcities = zeros(numCities, 2);
    individualcities(1,:) = cities(1,:);
    for j = 2:numCities
        individualcities(j,:) = cities(cityOrder(j - 1),:);
    end
    [individualPath, individualLength] = evaluateGraph(individualcities);
    P(i,:) = {individualPath, individualLength};
end
Np = round(N / 5) * 2; % Parent size
Pc = 0.8; % Crossover probability
Pm = 1 / Np; % Mutation probability

% Evaluate population
P = sortrows(P, 2);

% create the initial plot
figure
p = plot(path(:, 1), path(:, 2), 'o-');
% store gca so we can change the title later
ha = gca;
xlim([0, gridSize])
ylim([0, gridSize])
title(['Path length: ', num2str(length)])
% set the data sources for the plot so we can update the plot later
p.XDataSource = 'P{1,1}(:, 1)';
p.YDataSource = 'P{1,1}(:, 2)';

% Algorithm execution
while iterations <= maxIterations
    P_new = cell(N, 2); % The next generation

    % Crossover parents
    for i = 1:N
        currentParrent = P(mod(i, Np) + 1,:);
        r = rand;
        if r < Pc
            currentParrent2 = P(mod(i+(Np/2), Np) + 1,:);
            parrent1Cities = currentParrent{1}(1:numCities,:);
            parrent2Cities = currentParrent2{1}(1:numCities,:);
            cut = randi([round(numCities * 0.15), round(numCities * 0.75)]);
            childCities = parrent1Cities(1:(cut - 1),:);
            childCities = cat(1, childCities, parrent2Cities);
            childCities = unique(childCities, 'rows', 'stable');
            [childPath, childLength] = evaluateGraph(childCities);
            P_new(i,:) = {childPath, childLength};
        else
            P_new(i,:) = currentParrent;
        end
    end
    
    % Mutate children
    for i = 1:N
        r = rand;
        if r < Pm
            mutatedCities = randomSwap(P_new{i,1}(1:numCities,:));
            [mutatedPath, mutatedLenth] = evaluateGraph(mutatedCities);
            P_new(i,:) = {mutatedPath, mutatedLenth};
        end
    end

    % Adjust mutation rate
    numberOfDuplicates = width([P_new{:,2}]) - width(unique([P_new{:,2}]));
    Pm = numberOfDuplicates / N;

    % Evaluate population
    P = P_new;
    P = sortrows(P, 2);

    iterations = iterations + 1;
end

figure
subplot(1, 2, 1)
plot(initialPath(:, 1), initialPath(:, 2), 'o-')
xlim([0, gridSize])
ylim([0, gridSize])
title(['Initial random path length: ', num2str(initialLength)])

subplot(1, 2, 2)
plot(P{1,1}(:, 1), P{1,1}(:, 2), 'o-')
xlim([0, gridSize])
ylim([0, gridSize])
title(['Final path length: ', num2str(P{1,2})])

function newPath = randomSwap(oldPath)
    % copy the old path to the new path
    newPath = oldPath;

    % generate two random indices that are not the first
    a = randi([2, size(newPath, 1)]);
    b = randi([2, size(newPath, 1)]);
        
    % swap the cities at the two indices
    newPath([a, b], :) = newPath([b, a], :);
end