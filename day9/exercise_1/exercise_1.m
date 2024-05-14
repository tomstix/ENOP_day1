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

% Parameters
N = 100; % Population size
k_max = 300; % Maximum number of iterations
p_c = 0.8; % Crossover probability
p_m = 0.1; % Mutation probability
n_elite = 40; % Number of elite individuals
sigma_s = 0.04; % Sharing function threshold
alpha_s = 1; % Sharing function shape parameter
sigma_m = 0.1; % Crossover threshold
d_p_max = 1; % Maximum length of the pareto front
pause_time = 0.2;

rng default

% Initialization
X = lim(:,1) + (lim(:,2) - lim(:,1)) .* rand(dim, N);
F = fn(X);
D = zeros(1, N);
for i = 1:N
    D(i) = count_dominated_by(F(:,i), F);
end
SF = compute_shared_fitnesses(X, D, sigma_s, alpha_s);
X_elite = X(:,1:n_elite); % Might want to sort before selecting elite
% Initial Plot
figure;
pl = scatter(F(1,:), F(2,:), 'filled');
% xlim(plot_limits(1,:));
% ylim(plot_limits(2,:));

% Main loop
for k = 1:k_max
    % calculate sigma_share
    F_z = [min(F(1,:)), min(F(2,:))];
    % find the two points with the biggest distance in F
    combs = nchoosek(1:N, 2);
    max_dist = 0;
    best_comb = [0, 0];
    for i = 1:size(combs, 1)
        d = norm(F(:, combs(i,1)) - F(:, combs(i,2)));
        if d > max_dist
            max_dist = d;
            best_comb = combs(i,:);
        end
    end
    if best_comb == [0, 0]
        error("No best combination found");
    end
    F_x = F(:, best_comb(1));
    F_y = F(:, best_comb(2));
    if F_x(1) > F_y(1)
        temp = F_x;
        F_x = F_y;
        F_y = temp;
    end
    d1 = norm(F_x - F_z);
    d2 = norm(F_y - F_z);
    d_min = sqrt(d1^2 + d2^2);
    d_max = d1 + d2;
    d_k = (d_min + d_max) / 2;
    if d_k > d_p_max
        d_p_max = d_k;
    end
    sigma_share = N^(1/1-m) * d_k;
    sigma_mate = 3 * N^(1/1-m) * d_k;

    % Selection
    p_idx = tournament_selection(D, SF);
    X_new = crossover(X(:,p_idx), F(:, p_idx), p_c, sigma_m);
    X_new = mutation(X_new, p_m, lim);
    X = [X_new, X_elite];
    F = fn(X);
    D = zeros(1, size(X, 2));
    for i = 1:size(X, 2)
        D(i) = count_dominated_by(F(:,i), F);
    end
    SF = compute_shared_fitnesses(X, D, sigma_s, alpha_s);
    [~, idx] = sort(SF, "descend");
    X = X(:,idx(1:N));
    F = F(:,idx(1:N));
    D = D(idx(1:N));
    SF = SF(idx(1:N));
    X_elite = X(:,1:n_elite);
    % Plot
    hold on;
    delete(pl);
    if exist('pl2', 'var')
        delete(pl2);
    end
    pl = plot(F(1,:), F(2,:), 'bo');
    pl2 = plot(F(1,1:n_elite), F(2,1:n_elite), 'ro');
    xlim(plot_limits(1,:));
    ylim(plot_limits(2,:));
    title(['Iteration ', num2str(k)]);
    drawnow limitrate;

    fprintf("Iteration %d \t d_k: %f \t sigma_share: %f \t sigma_mate: %f\n", k, d_k, sigma_share, sigma_mate)

    pause(pause_time);
end

function d = count_dominated_by(f, set)
% Count the number of solutions in set that are dominating f
d = 0;
for i = 1:size(set, 2)
    if all(f >= set(:,i)) && any(f > set(:,i))
        d = d + 1;
    end
end
end

function SF = compute_shared_fitnesses(X, D, sigma, alpha)
% Compute the shared fitness of a solution
% X: population
% F: objective values of the population
% D: domination count of the population
% sigma: sharing threshold
% alpha: sharing function shape parameter
N = size(X, 2);
SF = zeros(1, N);
for i = 1:N
    SF(i) = 0;
    for j = 1:N
        d = norm(X(:,i) - X(:,j));
        if d < sigma
            SF(i) = SF(i) + (1 - d/sigma)^alpha;
        end
    end
end
SF = (1 ./ (1 + D)) ./ SF;
end

function idx = tournament_selection(D, SF)
% Perform tournament selection
% D: dominated by count of the candidates
% SF: shared fitness of the candidates
N = size(D, 2);
idx = zeros(1, N);
for i = 1:2*N
    r = randi(N, 1, 2);
    c1_dominated = D(r(1)) > 0;
    c2_dominated = D(r(2)) > 0;
    if c1_dominated == c2_dominated
        f1 = SF(r(1));
        f2 = SF(r(2));
        if f1 > f2
            idx(i) = r(1);
        else
            idx(i) = r(2);
        end
    elseif ~c1_dominated
        idx(i) = r(1);
    else
        idx(i) = r(2);
    end
end
end

function X_new = crossover(X, F, p_c, sigma_m)
% Perform crossover
% X: parents
% F: objective values of the parents
% p_c: crossover probability
% sigma_m: mating threshold
dim = size(X, 1);
N = size(X, 2);
X_new = zeros(dim, N/2);
for i = 1:N/2
    if rand < p_c
        p1 = i;
        p2 = 2 * i;
        d = norm(F(:, i) - F(:, 2*i));
        if d ~= 0 && d > sigma_m
            searchperm = randperm(N/2) + N/2;
            for j = searchperm
                d = norm(F(:, i) - F(:, j));
                if d <= sigma_m && d ~= 0
                    p2 = j;
                    break;
                end
            end
        end
        alpha = rand;
        X_new(:, i) = alpha * X(:, p1) + (1 - alpha) * X(:, p2);
    else
        % copy parent
        X_new(:, i) = X(:, i);
    end
end
end

function X_new = mutation(X, p_m, lim)
% Perform mutation
% X: population to mutate
% p_m: mutation probability
% lim: limits of the parameters
dim = size(X, 1);
N = size(X, 2);
X_new = X;
for i = 1:N
    if rand < p_m
        % perform mutation
        X_new(:, i) = lim(:,1) + (lim(:,2) - lim(:,1)) .* rand(dim, 1);
    end
end
end
