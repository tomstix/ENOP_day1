obstacle_probability = 0.3;
board_size = 20;

b = make_board(board_size, obstacle_probability);

% find all non-obstacle nodes at the initial row
start_nodes_x = find(b(1,:) ~= 1);

% matrix that stores whether a node has been explored
explored = zeros(board_size);

% add all start nodes to frontier
frontier = {};
for i = 1:size(start_nodes_x, 2)
    x = start_nodes_x(i);
    coords = {[x, 1]};
    frontier = [frontier, coords];
end

while numel(frontier) > 0
    node = frontier{1};
    explored(node(1), node(2)) = 1;
    if node(2) == board_size
        fprintf("Found goal node!")
        disp(node)
        break
    end
    frontier(1) = [];
    nb = get_neighbours(b, node(1), node(2));
    % add all unexplored neighbours to frontier
    for i = 1:size(nb, 2)
        if explored(nb{i}(1), nb{i}(2)) == 0
            frontier = [frontier, nb(i)];
        end
    end
end

function neighbours = get_neighbours(board, x, y)
    % create all possible moves
    board_size = size(board, 2);
    neighbours = {};
    if (x+1) > 1 && (x+1) <= board_size && board((x+1), y) == 0
        neighbours = [neighbours, [(x+1), y]];
    end
    if (x-1) > 1 && (x-1) <= board_size && board((x-1), y) == 0
        neighbours = [neighbours, [(x-1), y]];
    end
    if (y+1) > 1 && (y+1) <= board_size && board(x, (y+1)) == 0
        neighbours = [neighbours, [x, (y+1)]];
    end
    if (y-1) > 1 && (y-1) <= board_size && board(x, (y-1)) == 0
        neighbours = [neighbours, [x, (y-1)]];
    end
end

function plot_board(board)
    figure
    hold on
    [r, c] = find(board);
    scatter(c, r);
    yline(0);
    yline(height(board))
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