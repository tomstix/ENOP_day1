
a = -10;
b = 10;
X = randi([a, b], 1, 5);
Y = randi([a, b], 1, 5);
pgon = polyshape(X, Y); % Create the initial polygon

% Calculate the centroid of the polygon
[xc, yc] = centroid(pgon);
pgon=translate(pgon,-xc,-yc);

figure;
axis equal;
axis([a-5 b+5 a-5 b+5]); % axis limits
hold on;
plot(0,0,'ro'); % Mark the center (0,0) with a red circle
pgonPlot = plot(pgon); % initial polygon plot


while true
    pgon = rotate(pgon, -5, [0, 0]); % Rotate the polygon by 5 degrees clockwise around its centroid
    pgonPlot.Shape = pgon; % plot update
    drawnow; % Update the plot figure
    pause(0.1); % Pause to make visible
end