clear; clc; close all;

fn_select = 2;
switch fn_select
    case 1
        fn = @ft;
        lim = [0.5, 1.0; -2.0, 2.0; -2.0, 2.0];
        dim = 3;
        m = 2;
        plot_limits = [0.4, 1; 1, 5];
    case 2
        fn = @ff;
        lim = [-2 * ones(8,1), 2 * ones(8,1)];
        dim = 8;
        m = 2;
        plot_limits = [0, 1; 0, 1];
end

p_c = 0.9;
p_m = 0.1;
mutation_beta = 3;
sigma_share = 0.5;
alpha = 1;
sigma_mate = 3 * sigma_share;
N = 100;
n_elite = round(N * 0.4);
p_hist = zeros(1, N);
p_break = 0.01;

num_iterations = 200;

P = create_population(N, lim);
P = [P, evaluate_function(fn, P)];
[shared_fitnesses, q] = calculate_shared_fitnesses(P, alpha, sigma_share);
P = [P, shared_fitnesses, q];
P = sortrows(P, 3, 'descend');
nondom_idx = find(cell2mat(P(:,4)) == 0);
nondoms = P(nondom_idx, :);

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
xlim(plot_limits(1,:));
ylim(plot_limits(2,:));

k = 0;
while k < num_iterations
    k = k + 1;
    
    % calculate sigma_share
    F_z = [min(F(:,1)), min(F(:,2))];
    % find the two points with the biggest distance in F
    combs = nchoosek(1:N, 2);
    max_dist = 0;
    best_comb = [0, 0];
    for i = 1:size(combs, 1)
        d = norm(F(combs(i,1),:) - F(combs(i,2),:));
        if d > max_dist
            max_dist = d;
            best_comb = combs(i,:);
        end
    end
    F_x = F(best_comb(1),:);
    F_y = F(best_comb(2),:);
    d1 = norm(F_x - F_z);
    d2 = norm(F_y - F_z);
    d_min = sqrt(d1^2 + d2^2);
    d_max = d1 + d2;
    d_k = (d_min + d_max) / 2;
    sigma_share = N^(1/1-m) * d_k;
    sigma_mate = 3 * sigma_share;

    pll1 = plot([F_x(1), F_z(1)], [F_x(2), F_z(2)], 'g');
    pll2 = plot([F_y(1), F_z(1)], [F_y(2), F_z(2)], 'g');
    pll3 = plot([F_x(1), F_y(1)], [F_x(2), F_y(2)], 'g');
    
    parents = parent_selection(P);
    offspring = crossover(parents, p_c, sigma_mate);
    offspring = mutation(offspring, p_m, mutation_beta, lim);
    offspring = mat2cell(offspring, ones(1,N), height(lim));
    offspring = [offspring, evaluate_function(fn, offspring), cell(N,1), cell(N,1)];
    P = [elite; offspring];
    [shared_fitnesses, q] = calculate_shared_fitnesses(P, alpha, sigma_share);
    P(:,3) = shared_fitnesses;
    P(:,4) = q;
    P = sortrows(P, 3, 'descend');
    P = P(1:N, :);
    elite = P(1:n_elite, :);
    F = cell2mat(P(:,2));
    E = cell2mat(elite(:,2));
    title(['Iteration ', num2str(k)]);
    refreshdata
    drawnow limitrate
    fprintf('Iteration %d, sigma_s: %.4f, d_k = %0.3f\n', k, sigma_share, d_k);
    delete(pll1);
    delete(pll2);
    delete(pll3);

    if k < 10
        continue;
    end
    nondom_idx_new = find(cell2mat(P(:,4)) == 0);
    nondoms_new = P(nondom_idx_new, :);
    n_dash = 0;
    for i = 1:height(nondoms_new)
        for j = 1:height(nondoms)
            if is_dominated_by(nondoms{j,2}, nondoms_new{i,2})
                n_dash = n_dash + 1;
                break;
            end
        end
    end
    p = n_dash / length(nondom_idx_new);
    p_hist(k) = p;
    p_avg = mean(p_hist(k-9:k));
    fprintf('p: %.4f, p_avg: %4f\n', p, p_avg);
    if p_avg < p_break
        break;
    end
    nondoms = nondoms_new;

    fprintf("\n");
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

function SF = sharing_function(f1, f2, alpha, sigma_share)
d = norm(f1 - f2);
if d < sigma_share
    SF = 1 - (d / sigma_share)^alpha;
else
    SF = 0;
end
end

function SF = shared_fitness_sum(f, F, alpha, sigma_share)
N = height(F);
SF = 0;
for i = 1:N
    SF = SF + sharing_function(f, F(i,:), alpha, sigma_share);
end
end

function [f, q] = calculate_shared_fitnesses(P, alpha, sigma_share)
% P has to be a cell array with the individuals in the first row and function values in the second row
N = height(P);
F = cell2mat(P(:,2));
f = zeros(N, 1);
q = zeros(N, 1);
for i = 1:N
    q(i) = is_dominated(F(i,:), F);
    r = 1 + q(i);
    f(i) = 1 / r;
    sum_sf = shared_fitness_sum(F(i,:), F, alpha, sigma_share);
    f(i) = f(i) / sum_sf;
end
f = num2cell(f);
q = num2cell(q);
end

function S = parent_selection(P)
P_orig = P;
S = [P_orig; P_orig];
SF = cell2mat(P(:,3));
P = cell2mat(P(:,1));
D = cell2mat(P_orig(:,4));
N = height(P);
for i = 1:2*N
    c = randi(N, 1, 2);
    c1_dominated = D(c(1)) > 0;
    c2_dominated = D(c(2)) > 0;
    if c1_dominated == c2_dominated
        f1 = SF(c(1));
        f2 = SF(c(2));
        if f1 > f2
            S(i, :) = P_orig(c(1), :);
        else
            S(i, :) = P_orig(c(2), :);
        end
    else
        if c1_dominated == false
            S(i, :) = P_orig(c(1), :);
        else
            S(i, :) = P_orig(c(2), :);
        end
    end
end
end

function C = crossover(P, p_c, sigma_mate)
F = cell2mat(P(:,3));
P = cell2mat(P(:,1));
N = height(P);
dim = width(P);
C = zeros(N/2, dim);
for i = 1:N/2
    p1 = P(2*i-1,:);
    f1 = F(2*i-1,:);
    if rand() < p_c
        found = false;
        % search for mate within sigma_mate
        for j = 1:N
            if i ~= j
                p2 = P(j,:);
                f2 = F(j,:);
                d = norm(f1 - f2);
                if d < sigma_mate && d ~= 0
                    found = true;
                    break;
                end
            end
        end
        % if no mate was found, choose random mate
        if ~found
            p2 = P(ceil(N*rand()),:);
        end
        % intermediate crossover
        for k = 1:dim
            alpha = rand();
            C(i,k) = alpha*p1(k) + (1-alpha)*p2(k);
        end
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