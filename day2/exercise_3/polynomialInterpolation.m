function polynomialInterpolation(myMatrix)
   % First check if the argument is valid
   if size(myMatrix, 2) ~= 2
       warning("Invalid argument! Should be matrix of size nx2")
       return
   end

   % Define a few variables needed to calculate the polynomial
   n = size(myMatrix, 1);
   N = n - 1;
   Y = myMatrix(:,2);
   A = zeros(n, N + 1);
   for j = 1:n
       for i = 1:N + 1
           A(j,i) = myMatrix(j,1)^(N - (i - 1));
       end
   end

   % Calculate the coefficients of the polynomial
   B = (A' * A)^-1 * A' * Y;
   polynomial = B';

   % Find the derivative
   derivative = polyder(polynomial);

   % Find the maxima and minima
   derivativeRoots = roots(derivative);
   minMax = zeros(size(derivativeRoots, 1),2);
   for i = 1:size(minMax, 1)
       minMax(i, 1) = derivativeRoots(i);
       minMax(i, 2) = polyval(polynomial,derivativeRoots(i));
   end

   % Plot everything
   x = min(myMatrix(:,1)):.1:max(myMatrix(:,1));
   figure
   set(0, "DefaultLineLineWidth", 1);
   plot(x,polyval(polynomial,x));
   hold on;
   scatter(myMatrix(:,1), myMatrix(:,2), "LineWidth", 1);
   plot(x, 0*x);
   plot(x, polyval(derivative,x)) % Derivative
   scatter(minMax(:,1), minMax(:,2), 100, "x", "LineWidth", 1.5) % Minima and maxima
   legend('Interpolation polynomial {\it p}', 'Data points', 'Zero level', 'First-order derivative of {\it p}', 'Locations of local maxima/minima', "Location","southeast")
   xlim([min(myMatrix(:,1)) max(myMatrix(:,1))]);
   grid on;
   hold off;
end