numRectangles = 70; % the number of rectangles to pack
binWidth = 10; % the width of the bin
binHeight = 200; % the height of the bin, might need to be increased for more rectangles
rectangleMaxWidth = 8; % the maximum width of a rectangle
rectangleMaxHeight = 6; % the maximum height of a rectangle

% generate random rectangles
rectangles = struct('size', {}, 'position', {}, 'color', {});
for i = 1:numRectangles
    rectangles(i).size = [randi(rectangleMaxWidth-1)+1 randi(rectangleMaxHeight-1)+1];
    rectangles(i).position = [0 0];
    rectangles(i).color = rand(1,3);
end

% create new figure for the updating plot
figure;
ax = axes;

% this is to break the loop if we are stuck
numUnsuccessfulSwaps = 0;
maxUnsuccessfulSwaps = 10000;

% initial packing
[rectangles, maxHeight] = packRectangles(rectangles, binWidth, binHeight, false);
initialRectangles = rectangles;
initialMaxHeight = maxHeight;
drawRectangles(ax, rectangles, maxHeight, binWidth, binHeight);
% pause a bit to see the progress better
pause(0.5)

% repeat until we are stuck
while numUnsuccessfulSwaps < maxUnsuccessfulSwaps
    % swap two random rectangles in the list
    newRectangles = swapRectangles(rectangles);

    % pack the new rectangles and calculate the new max height
    [newRectangles, newMaxHeight] = packRectangles(newRectangles, binWidth, binHeight, true);
    
    % update the rectangles if the new packing is higher
    if newMaxHeight < maxHeight
        rectangles = newRectangles;
        maxHeight = newMaxHeight;
        numUnsuccessfulSwaps = 0; % reset the counter
        drawRectangles(ax, rectangles, maxHeight, binWidth, binHeight); % update the plot
    else
        numUnsuccessfulSwaps = numUnsuccessfulSwaps + 1; % increase the counter
    end
end

% plot final comparison
figure;
f = subplot(1, 2, 1);
drawRectangles(f, initialRectangles, initialMaxHeight, binWidth, binHeight)
title(sprintf('Initial packing\nMax height: %d', initialMaxHeight))

f = subplot(1, 2, 2);
drawRectangles(f, rectangles, maxHeight, binWidth, binHeight)
title(sprintf('Final packing\nMax height: %d', maxHeight))

% function to draw the rectangles / updates figure f
function drawRectangles(ax, rectangles, maxHeight, binWidth, binHeight)
cla(ax); % clear the plot
% axis equal % make the axes equal. This might look bad for a high height to width ratio
% iterate over the rectangles and draw them
for i = 1:size(rectangles, 2)
    xCoord = rectangles(i).position(1);
    yCoord = rectangles(i).position(2);
    width = rectangles(i).size(1);
    height = rectangles(i).size(2);
    rectangle('Position',[xCoord-1 yCoord-1 width height],'FaceColor',rectangles(i).color,'EdgeColor','b', 'LineWidth',1)
end
xlim([0 binWidth])
ylim([0 binHeight + 10])
title(['Max height: ', num2str(maxHeight)])
drawnow % update the plot immediately, don't wait for the script to finish
pause(0.5) % pause a bit to see the progress better
end

function [rectangles, maxHeight] = packRectangles(rectangles, binWidth, binHeight, autoRotate)
bin = zeros(binWidth,binHeight); % initialize the bin with zeros, so it is empty
maxHeight = 0;
for i = 1:size(rectangles, 2) % iterate over the rectangles and place them individually
    rect = rectangles(i);
    [bin, rect] = placeRectangle(bin, rect, autoRotate);
    rectangles(i) = rect;
    yCoord = rect.position(2);
    height = rect.size(2);
    if yCoord + height > maxHeight % update the max height if needed
        maxHeight = yCoord + height;
    end
end
end

function [bin, rect] = placeRectangle(bin, rect, autoRotate)
% iterate over the bin, bottom to top, left to right
for y = 1:width(bin)
    for x = 1:height(bin)
        if bin(x, y) == 0 % check if the current position is empty
            if x + rect.size(1) - 1 <= height(bin) % check if the rectangle fits in the bin
                if sum(sum(bin(x:x+rect.size(1)-1, y:y+rect.size(2)-1))) == 0 % check if the rectangle fits without overlapping
                    rect.position = [x, y];
                    bin(x:x+rect.size(1)-1, y:y+rect.size(2)-1) = 1;
                    return
                end
            elseif x + rect.size(2) - 1 <= height(bin) && autoRotate % check if the rectangle fits rotated
                if sum(sum(bin(x:x+rect.size(2)-1, y:y+rect.size(1)-1))) == 0 % check if the rectangle fits without overlapping
                    % swap the width and height of the rectangle
                    rect.size = [rect.size(2), rect.size(1)];
                    rect.position = [x, y];
                    bin(x:x+rect.size(1)-1, y:y+rect.size(2)-1) = 1;
                    return
                end
            end
        end
    end
end
end

function rectangles = swapRectangles(rectangles)
% generate two random indices
i = randi([1, size(rectangles, 2)]);
j = randi([1, size(rectangles, 2)]);

% swap the rectangles at the two indices
temp = rectangles(i);
rectangles(i) = rectangles(j);
rectangles(j) = temp;
end