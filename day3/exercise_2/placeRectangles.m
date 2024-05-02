function [rectangles, bin] = placeRectangles(rectangles, bin)
    % Rezero bin
    bin = zeros(height(bin),width(bin));
    
    % Create cursor
    cursor = [0 0];

    for i = 1:numel(rectangles)
        rectangles(i).position = {};
        while isempty(rectangles(i).position)
            [cursor(1), cursor(2)] = find(bin==0,1,'first');
            if cursor(1) - 1 + rectangles(i).size(1) <= height(bin) && sum(sum(bin(cursor(1):cursor(1)-1+rectangles(i).size(1), cursor(2):cursor(2)-1+rectangles(i).size(2)))) == 0
                rectangles(i).position = cursor - 1;
                bin(cursor(1):cursor(1)-1+rectangles(i).size(1), cursor(2):cursor(2)-1+rectangles(i).size(2)) = 1;
                bin(bin==2) = 0;
                break
            elseif cursor(1) - 1 + rectangles(i).size(2) <= height(bin) && sum(sum(bin(cursor(1):cursor(1)-1+rectangles(i).size(2), cursor(2):cursor(2)-1+rectangles(i).size(1)))) == 0
                rectangles(i).size = [rectangles(i).size(2) rectangles(i).size(1)];
                rectangles(i).position = cursor - 1;
                bin(cursor(1):cursor(1)-1+rectangles(i).size(1), cursor(2):cursor(2)-1+rectangles(i).size(2)) = 1;
                bin(bin==2) = 0;
                break
            else
                bin(cursor(1), cursor(2)) = 2;
            end
        end
    end
end