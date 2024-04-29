obstacle_probability = 0.3;
size = 20;

b = make_board(size, obstacle_probability);



function neighbours = get_neighbours(board, x, y)
    % create all possible moves
    neighbours{1} = [x+1, y];
    neighbours{2} = [x, y+1];
    neighbours{3} = [x-1, y];
    neighbours{4} = [x, y-1];
    % delete all invalid moves
    for i = 1:numel(neighbours)
        if neighbours{i}(1) > width(board) || ...
           neighbours{i}(2) > height(board) || ...
           neighbours{i}(1) < 0 || neighbours{i}(2) < 0 || ...
           board(neighbours{i}(1), neighbours{i}(2)) == 1
            neighbours{i} = [];
        end
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