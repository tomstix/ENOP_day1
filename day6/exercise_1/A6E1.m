%Exercise 1

results_unc = zeros(length(n), 2);  % [fval, funcCount]
results_search = zeros(length(n), 2); % [fval, funcCount]
n = [3, 5, 10];
%options = optimset('Display', 'iter', 'TolFun', 10^-10);
options = optimset('Display','iter','TolFun',10^-10,'OutputFcn',@optimplotfval,'MaxFunEvals',5000,'MaxIter',5000); % For case n=10 for fminsearch

for i = 1:length(n)
    % Rosenbrock function for each dimension n
    f = @(x) sum((1-x(1:end-1)).^2 + 100.*(x(2:end)-x(1:end-1).^2).^2);
    
    x0 = zeros(n(i), 1); % Start vector
    disp(['Case n = ', num2str(n(i))]);
    
    % Using fminunc
    [x_unc, fval_unc, exitflag_unc, output_unc] = fminunc(f, x0, options);
    disp('Results using fminunc:');
    disp(['Minimum found: ', num2str(fval_unc)]);
    disp(['Function evaluations: ', num2str(output_unc.funcCount)]);
    results_unc(i,:) = [fval_unc, output_unc.funcCount];  % Store results
    figure;

    % Using fminsearch
    [x_search, fval_search, exitflag_search, output_search] = fminsearch(f, x0, options);
    disp('Results using fminsearch:');
    disp(['Minimum found: ', num2str(fval_search)]);
    disp(['Function evaluations: ', num2str(output_search.funcCount)]);
    results_search(i,:) = [fval_search, output_search.funcCount]; % Store results
    figure;
end

disp("__________________Results Overview________________________")
disp("n:"+n(1)+"  fminunc:"+results_unc(1,1)+"   f_eval:"+results_unc(1,2)+"  fminsearch:"+results_search(1,1)+"  fun_eval:"+results_search(1,2)+"")
disp("n:"+n(2)+"  fminunc:"+results_unc(2,1)+"  f_eval:"+results_unc(2,2)+"  fminsearch:"+results_search(2,1)+"  fun_eval:"+results_search(2,2)+"")
disp("n:"+n(3)+" fminunc:"+results_unc(3,1)+"  f_eval:"+results_unc(3,2)+"  fminsearch:"+results_search(3,1)+"  fun_eval:"++results_search(3,2)+"")
