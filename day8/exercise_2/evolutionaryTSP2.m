% version 2 of the evolutionary algorithm for the TSP
function [history, final_P] = evolutionaryTSP2(cities, num_individuals, max_evals, p_c, p_m_0, print, plot, grid_size)
if nargin < 7
    print = false;
end
if nargin < 8
    plot = false;
end
num_cities = size(cities, 1);
p_m = p_m_0; % initial mutation probability
P = initializePopulation(cities, num_individuals);
for i = 1:num_individuals
    if width(unique(P{i, 1})) ~= num_cities % check for duplicate cities
        error('Duplicate cities in path');
    end
    P{i, 2} = evaluatePath(cities, P{i, 1}); % initial evaluation
end
P = sortrows(P, 2); % sort population by path length so we have the best path in the first row
best_path = P{1, 2};

history = zeros(max_evals/num_individuals, 3);

if plot
    figure
    plotPath(gca, cities, P{1, 1}, P{1, 2}, grid_size);
    drawnow
end

k = 0;
num_evals = num_individuals;
% Execute as long as we're within the computational budget
while num_evals < max_evals 
    k = k + 1;
    % choose parents
    parents = selectParents(P);
    
    % crossover
    % we don't crossover the best path
    for i = 2:num_individuals
        parent1 = parents{2*i-1};
        parent2 = parents{2*i};
        % crossover with probability p_c
        if rand < p_c
            P{i, 1} = crossover(parent1, parent2);
            P{i, 2} = evaluatePath(cities, P{i, 1});
            num_evals = num_evals + 1;
        end
    end
    
    % mutation
    % we don't mutate the best path
    for i = 2:num_individuals
        % mutate with probability p_m
        if rand < p_m
            P{i, 1} = randomSwap(P{i, 1});
            P{i, 2} = evaluatePath(cities, P{i, 1});
            num_evals = num_evals + 1;
        end
    end
    
    % sort again
    P = sortrows(P, 2);
    
    % evaluate the diversity of the population
    S = std(cell2mat(P(:, 2)));
    % below a certain threshold we increase the mutation probability, above we decrease it
    if S < 100
        p_m = min(p_m*2, 0.9);
    else
        p_m = max(p_m/2, 0.1);
    end
    
    % update best path
    if P{1, 2} < best_path
        best_path = P{1, 2};
        if plot
            plotPath(gca, cities, P{1, 1}, P{1, 2}, grid_size);
            drawnow limitrate;
        end
    end
    if print
        fprintf('Generation: %d, Evals: %d, Best path length: %f, p_m: %f, sigma: %f\n', k, num_evals, P{1, 2}, p_m, S);
    end
    % this is necessary because it might happen that we don't evaluate anything, which will mess up the interpolation for the experiments
    if k > 1 && num_evals == history(k-1, 2)
        num_evals = num_evals + 1;
    end
    history(k, :) = [k, num_evals, P{1, 2}];
end
history = history(1:k, :);
% order cities according to the best path
final_P = cities(P{1, 1}, :);
end

function P = initializePopulation(cities, num_individuals)
num_cities = size(cities, 1);
P = cell(num_individuals, 2);
for i = 1:num_individuals
    P{i,1} = randperm(num_cities);
    P{i,2} = evaluatePath(cities, P{i,1});
end
end

function length = evaluatePath(cities, path)
length = 0;
for i = 1:width(path) - 1
    length = length + norm(cities(path(i), :) - cities(path(i + 1), :));
end
% add the distance from the last city back to the first city
length = length + norm(cities(path(end), :) - cities(path(1), :));
end

function child = crossover(parent1, parent2)
length = width(parent1);
% generate random crossover point
crossover_point = randi([1, length - 1]);
p1slice = parent1(1:crossover_point);
p2slice = parent2(crossover_point+1:end);
child = [p1slice, p2slice];
all = 1:length;
% find missing cities
missing = setdiff(all, child);
% find indices of duplicate cities
[~, w] = unique( child, 'stable' );
duplicate_indices = setdiff( 1:numel(child), w );
% replace duplicate cities with missing cities
child(duplicate_indices) = missing;
if width(unique(child)) ~= length
    error('Crossover failed');
end
end

function newPath = randomSwap(oldPath)
% copy the old path to the new path
newPath = oldPath;

% generate two random indices
i = randi([1, size(newPath, 2)]);
j = randi([1, size(newPath, 2)]);

% swap the cities at the two indices
newPath([i, j]) = newPath([j, i]);
end

function R = selectParents(P)
% this works as a tournament selection as described in the lecture
N = length(P);
R = cell(N, 1);
for i = 1:2*N
    c = ceil(N*rand(1,2));
    if P{c(1), 2} < P{c(2), 2}
        R{i} = P{c(1), 1};
    else
        R{i} = P{c(2), 1};
    end
end
end
