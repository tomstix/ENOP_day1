function binPacking(numRectangles, binWidth, binHeight, rectangleMaxWidth, rectangleMaxHeight, optimizationIterations)
    % generate random rectangles
    rectangles = struct('size', {}, 'position', {}, 'color', {});
    for i = 1:numRectangles
        rectangles(i).size = [randi(rectangleMaxWidth-1)+1 randi(rectangleMaxHeight-1)+1];
        rectangles(i).color = rand(1,3);
    end
    
    % Create bin
    bin = zeros(binWidth,binHeight);
    
    % Place rectangles
    [rectangles, bin] = placeRectangles(rectangles, bin);
    
    % Draw rectangles
    figure
    ax = axes;
    drawRectangles(rectangles, bin, ax);
    
    % Evaluate bin
    [x, y] = find(bin==1,1,'last');
    stackHight = y;
    
    % Optimize bin
    for i = 1:optimizationIterations
        % Replace rectangles
        switchPair = [randi(numRectangles) randi(numRectangles)];
        rectangles([switchPair(1) switchPair(2)]) = rectangles([switchPair(2) switchPair(1)]);
        [rectangles, bin] = placeRectangles(rectangles, bin);
    
        % Evaluate bin
        [x, y] = find(bin==1,1,'last');
        newStackHight = y;
        if newStackHight < stackHight
            stackHight = newStackHight;
            drawRectangles(rectangles, bin, ax);
        else
            rectangles([switchPair(1) switchPair(2)]) = rectangles([switchPair(2) switchPair(1)]);
        end
    end
end