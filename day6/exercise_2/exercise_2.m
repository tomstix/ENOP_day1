clear; clc;

n = 2; % change this to the wanted test case

x0 = zeros(1,n); % starting point

options = optimset('Display','iter'); % Set tolerances and display output
[x,f] = fmincon(@fn,x0,[],[],[],[],[],[],@normcon, options);
disp(['final result: [',num2str(x),']',]); disp(['final function value: ',num2str(f)]);

% show a plot in the 2D case
if n == 2
    fsurf(@(x,y) fn([x y]), [-1 1 -1 1]);
    hold on
    plot3(x(1),x(2),f,'r*');
    % sphere
    % axis equal
end

% constraint function for ||x|| <= 1
function [c, ceq] = normcon(x)
c = norm(x)^2 - 1;
ceq = [];
end

% actual funtion
function f = fn(x)
inner = 0;
for i = 1:width(x)
    inner = inner + x(i)/i;
end
f = exp(inner);
end