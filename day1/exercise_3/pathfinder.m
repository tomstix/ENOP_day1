clear; clc;

% probabilities = 0.25:0.05:0.75;
% num_iterations = 20;
% board_size = 20;
% 
% result = [];
% for i = 1:size(probabilities, 2)
%     sum_positives = 0;
%     for j = 1:num_iterations
%         [goal_found, b, goal_path] = find_path(probabilities(i), board_size);
%         if goal_found
%             sum_positives = sum_positives + 1;
%         end
%     end
%     p = sum_positives / num_iterations;
%     result = [result; probabilities(i), p];
% end

[goal_found, b, goal_path] = find_path(0.3, 20);
board_size = width(b);
% Plot
figure
grid on
xticks(1:1:board_size)
yticks(1:1:board_size)
hold on
[r, c] = find(b);
scatter(r, c, 'filled');
yline(1);
yline(height(b))
if goal_found
    scatter(goal_path(:,1), goal_path(:,2), 60)
end
xlim([0 board_size+1])
ylim([0 board_size+1])

function [goal_found, board, goal_path] = find_path(probability, board_size)
    board = make_board(board_size, probability);
    
    % find all non-obstacle nodes at the initial row
    start_nodes_x = find(board(:,1) ~= 1);
    
    % matrix that stores whether a node has been explored
    explored = zeros(board_size);
    
    % add all start nodes to frontier
    frontier = {};
    for i = 1:size(start_nodes_x, 1)
        x = start_nodes_x(i);
        coords = {[x, 1]};
        % add dummy first line, delete later
        path = {[x,1]};
        new_node = [coords, path];
        frontier = [frontier, {new_node}];
        explored(x, 1) = 1;
    end
    
    % Do Breath-First Search
    while numel(frontier) > 0
        node = frontier{1}{1};
        path = frontier{1}{2};
        explored(node(1), node(2)) = 1;
        frontier(1) = [];
        if node(2) == board_size
            goal_path = [path; node];
            len = height(goal_path);
            fprintf("Found goal node with %d steps!\n", len)
            break
        end
        nb = get_neighbours(board, node(1), node(2));
        % add all unexplored neighbours to frontier
        for i = 1:size(nb, 2)
            if explored(nb{i}(1), nb{i}(2)) == 0
                new_path = [path; node];
                new_node = [nb(i), {new_path}];
                frontier = [frontier, {new_node}];
            end
        end
    end
    
    goal_found = true;
    if ~exist("goal_path", "var") || height(goal_path) < 2
        warning("Could not find a goal")
        goal_found = false;
        goal_path = [];
    end
end

function neighbours = get_neighbours(board, x, y)
    % create all possible moves
    board_size = size(board, 2);
    neighbours = {};
    if (y+1) > 1 && (y+1) <= board_size && board(x, (y+1)) == 0
        neighbours = [neighbours, [x, (y+1)]];
    end
    if (x+1) > 1 && (x+1) <= board_size && board((x+1), y) == 0
        neighbours = [neighbours, [(x+1), y]];
    end
    if (x-1) > 1 && (x-1) <= board_size && board((x-1), y) == 0
        neighbours = [neighbours, [(x-1), y]];
    end
    if (y-1) > 1 && (y-1) <= board_size && board(x, (y-1)) == 0
        neighbours = [neighbours, [x, (y-1)]];
    end
end

function board = make_board(size, prob)
    board = zeros(size);
    for i = 1:height(board)
        for j = 1:width(board)
            x=rand;
            if x<prob
              board(i,j) = 1;
            end
        end
    end
end