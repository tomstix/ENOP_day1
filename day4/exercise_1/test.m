clear; clc;

% Select the function to test
fn_select = 1;

eps = 1e-3;

fn{1} = @(x) (x-2)^2;
x{1} = [0, 10];

fn{2} = @(x) x^2 + 3*exp(-2*x);
x{2} = [0, 10];

fn = fn{fn_select};
x1 = x{fn_select}(1);
x2 = x{fn_select}(2);

[xmin, numIterations] = goldenRatioSearch(fn, x1, x2, eps);
