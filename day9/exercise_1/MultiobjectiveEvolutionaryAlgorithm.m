clear; clc; close all;

% Variable that decides which function is used
functionSelect = 1;
% Valid values are 1 and 2, which function is which can be seen below

% Algorithm setup
generation = 1;
maxGeneration = 100;
numIndividuals = 100; % Number of individuals for the evolutionary algorithms
P = cell(numIndividuals,4); % Population
Pc = 0.7; % Crossover probability for the evolutionary algorithms
Pm = 0.1; % Mutation probability for the evolutionary algorithms
myAlpha = 1;
myBeta = 3;
myShareSigma = 0.25; % Removing dynamic sharing looks better
myMateSigma = 3 * myShareSigma;
E = cell(maxGeneration * numIndividuals/2,4); % Elite

% Defining the functions
switch functionSelect
    case 1
        funck = @(x) [1 - exp(-sum((x(1:end) - 1/sqrt(8)).^2)), 1 - exp(-sum((x(1:end) + 1/sqrt(8)).^2))];
        n = 8;
        range = [-2 2];
        for i = 1:numIndividuals
            P{i,1} = (range(2) - range(1)).*rand(1,n) + range(1);
            P{i,2} = funck(P{i,1});
            P{i,3} = 1;
            P{i,4} = 1;
        end
    case 2
        funck = @(x) [x(1), (1/x(1)) * (1 + (x(2)^2 + x(3)^2)^0.25 * (sin(50 * (x(2)^2 + x(3)^2)^0.1)^2 + 1))];
        n = 3;
        range = [0.5 1; -2 2];
        for i = 1:numIndividuals
            P{i,1} = zeros(1,n);
            P{i,1}(1) = (range(1,2) - range(1,1)).*rand(1,1) + range(1,1);
            P{i,1}(2:3) = (range(2,2) - range(2,1)).*rand(1,2) + range(2,1);
            P{i,2} = funck(P{i,1});
            P{i,3} = 1;
            P{i,4} = 1;
        end
end

% Evaluate pareto ranking
P = paretoRanking(P);
P = fitnessGrading(P, myAlpha, myShareSigma);

% Sort population
P = sortrows(P, 4, 'descend');
for i = 1:size(E,1)
    E(i,:) = P(1,:);
end

figure
scatter(cellfun(@(x) x(1), P(:,2)), cellfun(@(x) x(2), P(:,2)));
hold on;
scatter(cellfun(@(x) x(1), E(:,2)), cellfun(@(x) x(2), E(:,2)), 'MarkerEdgeColor', [1 0 0]);
xlim([0,1]);
ylim([0,1]);
grid on;
drawnow;
fprintf('Generation: %d, Best fitness: %f\n', generation, P{1,4});

% Algorithm execution
while generation < maxGeneration
    P_new = cell(numIndividuals,4); % The next generation
    
    % choose parents
    parents = selectParents(P);

    % Crossover parents
    for i = 1:numIndividuals
        parent1 = parents{2*i-1,1};
        parent2 = parents{2*i,1};

        restOfPopulation = 1:size(parents,1);
        restOfPopulation = [restOfPopulation(2*i:end), restOfPopulation(1:2*i-1)]; 
        for j = restOfPopulation
            distance = sqrt((parents{j,2}(1)-parent1(1))^2 + (parents{j,2}(2)-parent1(2))^2);
            if distance < myMateSigma && distance ~= 0
                parent2 = parents{j,1};
            end
        end

        split = rand;
        r = rand;
        if r < Pc
            P_new{i,1} = split .* parent1 + (1 - split) .* parent2;
        else
            P_new{i,1} = parent1;
        end
    end
    
    % Mutate children
    for i = 1:numIndividuals
        r = rand;
        if r < Pm
            for j = 1:n
                r = rand();
                if r > 0.5
                    if functionSelect == 2
                        if j == 1
                            delta_x = (range(1,2) - P{i,1}(j)) * (2 * (r - 0.5))^myBeta;
                        else
                            delta_x = (range(2,2) - P{i,1}(j)) * (2 * (r - 0.5))^myBeta;
                        end
                    else
                        delta_x = (range(2) - P{i,1}(j)) * (2 * (r - 0.5))^myBeta;
                    end
                else
                    if functionSelect == 2
                        if j == 1
                            delta_x = (range(1,1) - P{i,1}(j)) * (2 * (0.5 - r))^myBeta;
                        else
                            delta_x = (range(2,1) - P{i,1}(j)) * (2 * (0.5 - r))^myBeta;
                        end
                    else
                        delta_x = (range(1) - P{i,1}(j)) * (2 * (0.5 - r))^myBeta;
                    end
                end
                P{i,1}(j) = P{i,1}(j) + delta_x;
            end
        end
    end
    
    % Adjust sigma
    pointsX = cellfun(@(x) x(1), P(:,2));
    pointsY = cellfun(@(x) x(2), P(:,2));
    points = [pointsX, pointsY];
    points = sortrows(points, 2, 'descend');
    point1 = points(1,:);
    points = sortrows(points, 1, 'descend');
    point2 = points(1,:);
    d1 = point1(2) - point2(2);
    d2 = point2(1) - point1(1);
    dmin = sqrt(d1^2 + d2^2);
    dmax = d1 + d2;
    d = (dmax + dmin) / 2;
    myShareSigma = numIndividuals^(1/1-2) * (d / 2); % Adding " * 50" made it look better
    myMateSigma = 3 * myShareSigma;

    % Evaluate pareto ranking
    for i = 1:numIndividuals
        P_new{i,2} = funck(P_new{i,1});
    end
    P_new = paretoRanking(P_new);
    P_new = fitnessGrading(P_new, myAlpha, myShareSigma);

    % Sort population
    P_new = sortrows(P_new, 4, 'descend');
    P = P_new;
    
    E(size(E,1)/maxGeneration + 1:end) = E(1:end - size(E,1)/maxGeneration);
    for i = 1:size(E,1)/maxGeneration
        E(i,:) = P(1,:);
    end
    
    generation = generation + 1;
    
    cla(gca)
    scatter(cellfun(@(x) x(1), P(:,2)), cellfun(@(x) x(2), P(:,2)));
    hold on;
    scatter(cellfun(@(x) x(1), E(:,2)), cellfun(@(x) x(2), E(:,2)), 'MarkerEdgeColor', [1 0 0]);
    xlim([0,1]);
    ylim([0,1]);
    grid on;
    drawnow limitrate;

    fprintf('Generation: %d, Best fitness: %f\n', generation, P{1,4});
end

function rankedPopulation = paretoRanking(unrankedPopulation)
    populationSize = size(unrankedPopulation, 1);
    rankedPopulation = unrankedPopulation;

    for i = 1:populationSize
        rankedPopulation{i,3} = 1;
        restOfPopulation = 1:populationSize;
        restOfPopulation(i) = []; 
        for j = restOfPopulation
            if rankedPopulation{i,2}(1) - rankedPopulation{j,2}(1) > 0 && rankedPopulation{i,2}(2) - rankedPopulation{j,2}(2) > 0
                rankedPopulation{i,3} = rankedPopulation{i,3} + 1;
            end
        end 
    end
end

function gradedPopulation = fitnessGrading(nonGradedPopulation, fAlpha, fSigma)
    populationSize = size(nonGradedPopulation, 1);
    gradedPopulation = nonGradedPopulation;

    for i = 1:populationSize
        gradedPopulation{i,4} = 1/gradedPopulation{i,3};
        
        sharing = 1;
        pointi = nonGradedPopulation{i,2};
        restOfPopulation = 1:populationSize;
        restOfPopulation(i) = []; 
        for j = restOfPopulation
            pointj = nonGradedPopulation{j,2};
            distance = sqrt((pointj(1)-pointi(1))^2 + (pointj(2)-pointi(2))^2);
            if distance < fSigma
                sharing = sharing + (1 - (distance / fSigma)^fAlpha);
            end
        end

        gradedPopulation{i,4} = gradedPopulation{i,4} / sharing;
    end
end

function R = selectParents(P)
    % this works as a tournament selection as described in the lecture
    N = size(P, 1);
    n = size(P, 2);
    R = cell(N, n);
    for i = 1:2*N
        c = ceil(N*rand(1,2));
        if P{c(1),4} > P{c(2),4}
            R(i,:) = P(c(1),:);
        else
            R(i,:) = P(c(2),:);
        end
    end
end