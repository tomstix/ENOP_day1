function [path, length] = evaluateGraph(cities)
    % add the first city to the end of the list
    path = [cities; cities(1, :)];

    % calculate the distance between each city and add it to the total length
    length = 0;
    for i = 1:size(path, 1) - 1
        length = length + sqrt(sum((path(i, :) - path(i + 1, :)).^2));
    end
end