clear; clc; close all;

fn_select = 2;
switch fn_select
    case 1
        fn = @ft;
        lim = [0.5, 1.0; -2.0, 2.0; -2.0, 2.0];
    case 2
        fn = @ff;
        lim = [-2 * ones(8,1), 2 * ones(8,1)];
end

p_c = 0.7;
p_m = 0.1;
mutation_beta = 1;
sigma_share = 0.1;
alpha = 1;
sigma_mate = 3 * sigma_share;
N = 100;
n_elite = N/5;

num_iterations = 200;

P = create_population(N, lim);
P = [P, evaluate_function(fn, P)];
shared_fitnesses = calculate_shared_fitnesses(P, alpha, sigma_share);
P = [P, shared_fitnesses];
P = sortrows(P, 3, 'descend');

elite = P(1:n_elite, :);

figure
F = cell2mat(P(:,2));
E = cell2mat(elite(:,2));
hold on
pl1 = plot(F(:,1), F(:,2), 'ob');
pl2 = plot(E(:,1), E(:,2), 'or');
pl1.XDataSource = 'F(:,1)';
pl1.YDataSource = 'F(:,2)';
pl2.XDataSource = 'E(:,1)';
pl2.YDataSource = 'E(:,2)';
xlim([0, 1]);
ylim([0, 1]);

k = 0;
while k < num_iterations
    k = k + 1;

    % calculate sigma_share
    d1 = max(F(:,1)) - min(F(:,1));
    d2 = max(F(:,2)) - min(F(:,2));
    d_max = d1 + d2;
    d_min = sqrt(d1^2 + d2^2);
    d_k = d_max - d_min/2;
    sigma_share = N^(1/1-0.5) * d_k/2;
    sigma_mate = 3 * sigma_share;

    parents = parent_selection(P);
    offspring = crossover(parents, p_c, sigma_mate);
    offspring = mutation(offspring, p_m, mutation_beta, lim);
    offspring = mat2cell(offspring, ones(1,N), height(lim));
    offspring = [offspring, evaluate_function(fn, offspring), cell(N,1)];
    P = [elite; offspring];
    shared_fitnesses = calculate_shared_fitnesses(P, alpha, sigma_share);
    P(:,3) = shared_fitnesses;
    P = sortrows(P, 3, 'descend');
    P = P(1:N, :);
    elite = P(1:n_elite, :);
    F = cell2mat(P(:,2));
    E = cell2mat(elite(:,2));
    % scatter(E(:,1), E(:,2), 50, 'r', 'filled');
    title(['Iteration ', num2str(k)]);
    refreshdata
    drawnow limitrate
    fprintf('Iteration %d, \\sigma_s: %.2f\n', k, sigma_share);
end

function f = evaluate_function(fn, x)
if iscell(x)
    f = fn(cell2mat(x)');
else
    f = fn(x');
end
f = num2cell(f', 2);
end

function d = is_dominated(f1, F)
d = 0;
for i = 1:height(F)
    if is_dominated_by(f1, F(i, :))
        d = d + 1;
    end
end
end

function d = is_dominated_by(f1, f2)
% true if f2 dominates f1
d = all(f2 <= f1) && any(f2 < f1);
end

function P = create_population(N, lim)
dim = height(lim);
P = zeros(N, dim);
for i = 1:dim
    P(:,i) = lim(i,1) + (lim(i,2) - lim(i,1)) * rand(1, N);
end
P = mat2cell(P, ones(1,N), dim);
end

function f = fitness(f, F)
q = is_dominated(f, F);
r = 1 + q;
f = 1 / r;
end

function SF = sharing_function(f1, f2, alpha, sigma_share)
d = norm(f1 - f2);
if d < sigma_share
    SF = 1 - (d/sigma_share)^alpha;
else
    SF = 0;
end
end

function f = shared_fitness(f, F, alpha, sigma_share)
N = height(F);
SF = 0;
for i = 1:N
    SF = SF + sharing_function(f, F(i,:), alpha, sigma_share);
end
f = fitness(f, F) / SF;
end

function f = calculate_shared_fitnesses(P, alpha, sigma_share)
% P has to be a cell array with the individuals in the first row and function values in the second row
N = height(P);
F = cell2mat(P(:,2));
f = zeros(N, 1);
for i = 1:N
    f(i) = shared_fitness(F(i,:), F, alpha, sigma_share);
end
f = num2cell(f);
end

function S = parent_selection(P)
SF = cell2mat(P(:,3));
F = cell2mat(P(:,2));
P = cell2mat(P(:,1));
N = height(P);
dim = width(P);
S = zeros(N, dim);
for i = 1:2*N
    c = ceil(N * rand(1, 2));
    c1 = P(c(1), :);
    c1_dominated = is_dominated(F(c(1), :), F) > 0;
    c2 = P(c(2), :);
    c2_dominated = is_dominated(F(c(2), :), F) > 0;
    if c1_dominated == c2_dominated
        f1 = SF(c(1));
        f2 = SF(c(2));
        if f1 > f2
            S(i, :) = c1;
        else
            S(i, :) = c2;
        end
    else
        if c1_dominated == false
            S(i, :) = c1;
        else
            S(i, :) = c2;
        end
    end
end
end

function C = crossover(P, p_c, sigma_mate)
N = height(P);
dim = width(P);
C = zeros(N/2, dim);
for i = 1:N/2
    p1 = P(i,:);
    if rand() < p_c
        found = false;
        % search for mate within sigma_mate
        for j = 1:N
            if i ~= j
                p2 = P(j,:);
                d = norm(p1 - p2);
                if d < sigma_mate
                    found = true;
                    break;
                end
            end
        end
        if ~found
            p2 = P(ceil(N*rand()),:);
        end
        % intermediate crossover
        alpha = rand();
        C(i,:) = alpha*p1 + (1-alpha)*p2;
    else
        % copy one parent
        C(i,:) = p1;
    end
end
end

function R = mutation(P, p_m, beta, lim)
N = height(P);
dim = width(P);
for j = 1:N
    if rand() < p_m
        for k = 1:dim
            r = rand();
            if r > 0.5
                delta_x = (lim(k,2) - P(j,k))*(2*(r - 0.5))^beta;
            else
                delta_x = (lim(k,1) - P(j,k))*(2*(0.5-r))^beta;
            end
            P(j,k) = P(j,k) + delta_x;
        end
    end
end
R = P;
end