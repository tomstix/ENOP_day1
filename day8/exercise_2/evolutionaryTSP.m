function [history, bestPath] = evolutionaryTSP(cities, maxEvaluations, N, Pc, Pm, print, plot, grid_size)
if nargin < 6
    print = false;
end
if nargin < 7
    plot = false;
end

% Algorithm setup
numCities = size(cities, 1);
evaluations = 1;
generation = 1;
history = zeros(maxEvaluations / 50, 3);
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
    evaluations = evaluations + 1;
    P(i,:) = {individualPath, individualLength};
end
Np = round(N / 5) * 2; % Parent size

% Sort population
P = sortrows(P, 2);

history(generation,:) = [generation, 1, P{1,2}];

if plot
    figure
    plotPath_alt(gca, cities, P{1, 1}, P{1, 2}, grid_size);
    drawnow
end

% Algorithm execution
while evaluations <= maxEvaluations
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
            evaluations = evaluations + 1;
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
            evaluations = evaluations + 1;
            P_new(i,:) = {mutatedPath, mutatedLenth};
        end
    end
    
    % Adjust mutation rate
    numberOfDuplicates = width([P_new{:,2}]) - width(unique([P_new{:,2}]));
    Pm = numberOfDuplicates / N;
    
    % Sort population
    P = P_new;
    P = sortrows(P, 2);
    
    if plot
        plotPath_alt(gca, cities, P{1, 1}, P{1, 2}, grid_size);
        drawnow limitrate;
    end
    
    generation = generation + 1;
    if evaluations == history(generation-1, 2)
        evaluations = evaluations + 1;
    end
    history(generation,:) = [generation, evaluations, P{1,2}];
    if print
        fprintf('Generation: %d, Evaluations: %d, Best length: %f\n', generation, evaluations, P{1,2});
    end
end
history = history(1:generation,:);
bestPath = P{1,1};
end

function newPath = randomSwap(oldPath)
% copy the old path to the new path
newPath = oldPath;

% generate two random indices that are not the first
a = randi([2, size(newPath, 1)]);
b = randi([2, size(newPath, 1)]);

% swap the cities at the two indices
newPath([a, b], :) = newPath([b, a], :);
end