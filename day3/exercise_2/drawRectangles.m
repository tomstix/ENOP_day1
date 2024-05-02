function drawRectangles(rectangles, bin, ax)
    cla(ax)
    %axis([0 height(bin) 0 width(bin)])
    for i = 1:numel(rectangles)
        rectangle('Position',[rectangles(i).position(1) rectangles(i).position(2) rectangles(i).size(1) rectangles(i).size(2)],'FaceColor',rectangles(i).color,'EdgeColor','b', 'LineWidth',1)
    end
    drawnow
end