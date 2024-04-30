[x, y] = meshgrid(linspace(-2, 2, 100), linspace(-2, 2, 100));

z=x.^2+y.^2;   % Equation for a circular paraboloid

% Points for tangency
x0=0.5;        
y0=-0.5;
z0=x0^2+y0^2;

z_tplane=-z0+1*(x-x0)-1*(y-y0); %Surface tangent plane

figure;
surf(x, y, z, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;

% Plotting the tangent plane
surf(x, y, z_tplane, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', 'r');

% Highlighting the point of tangency
plot3(x0, y0, z0, 'ko', 'MarkerSize', 5, 'LineWidth', 2);

% Labels and title
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Circular Paraboloid and its Tangent Plane at (0.5, -0.5)');
legend('Paraboloid', 'Tangent Plane', 'Point of Tangency');
