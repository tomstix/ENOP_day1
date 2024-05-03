function [xmin, numIterations] = goldenRatioSearch(f, x1, x2, eps)
t = (1+sqrt(5))/2;
numIterations = 0;
while true
    numIterations = numIterations + 1;
    
    % find x3 and x4 according to formula
    x3 = x2 - (x2 - x1)/t;
    x4 = x1 + (x2 - x1)/t;
    
    % if f(x3) > f(x4), the minimum is between x3 and x2, so we eliminate x1,
    % otherwise, the minimum is between x1 and x4, so we eliminate x2
    if f(x3) > f(x4)
        x1 = x3;
    else
        x2 = x4;
    end
    
    if abs(x2 - x1) < eps
        xmin = (x1 + x2)/2;
        break;
    end
end
end