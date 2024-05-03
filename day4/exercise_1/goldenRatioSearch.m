% Write a Matlab function that implements a golden ratio search algorithm.
% Let f be a function of a single parameter x.
% Let the (unique) minimum of f be bracketed by two points x1 and x2 (x1 < x2).
% Let t = (1+5^1/2)/2 (golden ratio).
% Find two points x3 and x4 so that (x2 – x3) = (x2 – x1)/t and (x4 – x1) = (x2 – x1)/t.
% Eliminate the point that has the largest value of f(xj), j = 1, 2, 3, 4.
% Set x1 = x3 (if x1 is eliminated) or x2 = x4 (if x2 is eliminated).
% Iterate the procedure until the minimum is located within the prescribed toletance.

function [xmin, numIterations] = goldenRatioSearch(f, x1, x2, eps)
fprintf("Searching for minimum of function %s\n", func2str(f))
t = (1+sqrt(5))/2;
numIterations = 0;
while true
    numIterations = numIterations + 1;
    fprintf('Cycle %d: x1 = %f, x2 = %f\n', numIterations, x1, x2)
    
    % find x3 and x4 according to formula
    x3 = x2 - (x2 - x1)/t;
    x4 = x1 + (x2 - x1)/t;
    
    % if f(x3) > f(x4), the minimum is between x3 and x2, so we eliminate x1,
    % otherwise, the minimum is between x1 and x4, so we eliminate x2
    if f(x3) > f(x4)
        fprintf('Changing x1 to %f\n', x3)
        x1 = x3;
    else
        fprintf('Changing x2 to %f\n', x4)
        x2 = x4;
    end

    fprintf('f(x1) = %f, f(x2) = %f\n', f(x1), f(x2))
    
    if abs(x2 - x1) < eps
        xmin = (x1 + x2)/2;
        fprintf('Minimum found at x = %f\n', xmin)
        break;
    end
    fprintf('\n')
end
end