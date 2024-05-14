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
k_max = 200; % Maximum number of iterations
p_c = 0.8; % Crossover probability
p_m = 0.3; % Mutation probability
n_elite = 40; % Number of elite individuals
alpha_s = 1; % Sharing function shape parameter
d_p_max = 1; % Maximum length of the pareto front
rand_beta = 2; % Beta for random mutation
pause_time = 0;

rng default

% Initialization
X = lim(:,1) + (lim(:,2) - lim(:,1)) .* rand(dim, N);
F = fn(X);
D = zeros(1, N);
for i = 1:N
    D(i) = count_dominated_by(F(:,i), F);
end
SF = compute_shared_fitnesses(X, D, 1, alpha_s);
X_elite = X(:,1:n_elite); % Might want to sort before selecting elite
% Initial Plot
figure;
hold on;
pl = scatter(F(1,:), F(2,:), 'filled', 'b');
pl_e = scatter(F(1,1:n_elite), F(2,1:n_elite), 'filled', 'r');
pl.XDataSource = 'F(1,:)';
pl.YDataSource = 'F(2,:)';
pl_e.XDataSource = 'F(1,1:n_elite)';
pl_e.YDataSource = 'F(2,1:n_elite)';
xlim(plot_limits(1,:));
ylim(plot_limits(2,:));

% Main loop
for k = 1:k_max
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
    d1 = abs(F_x(2) - F_y(2));
    d2 = abs(F_x(1) - F_y(1));
    d_min = sqrt(d1^2 + d2^2);
    d_max = d1 + d2;
    d_k = (d_min + d_max) / 2;
    if d_k > d_p_max
        d_p_max = d_k;
    end
    sigma_share = N^(1/1-m) * d_k;
    sigma_mate = 3 * N^(1/1-m) * d_k;
    
    % Selection
    X_new = crossover(X, F, D, SF, p_c, sigma_mate);
    X_new = mutation(X_new, p_m, lim, rand_beta);
    X = [X_new, X_elite];
    F = fn(X);
    D = zeros(1, size(X, 2));
    for i = 1:size(X, 2)
        D(i) = count_dominated_by(F(:,i), F);
    end
    SF = compute_shared_fitnesses(X, D, sigma_share, alpha_s);
    [~, idx] = sort(SF, "descend");
    X = X(:,idx(1:N));
    F = F(:,idx(1:N));
    D = D(idx(1:N));
    SF = SF(idx(1:N));
    X_elite = X(:,1:n_elite);
    % Plot
    refreshdata
    title(['Iteration ', num2str(k)]);
    drawnow limitrate;
    
    fprintf("Iteration %d \t d_k: %f \t sigma_share: %f \t sigma_mate: %f\n", k, d_k, sigma_share, sigma_mate)
    
    pause(pause_time);

    % rand_beta = rand_beta + rand_beta * (k / k_max);
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

function X_new = crossover(X, F, D, SF, p_c, sigma_mate)
dim = size(X, 1);
N = size(X, 2);
X_new = zeros(dim, N);
for i = 1:N
    % tournament selection for parent 1
    % select 2 parents randomly
    ps = randperm(N, 2);
    p1 = ps(1);
    p2 = ps(2);
    p1_dom = D(p1) > 0;
    p2_dom = D(p2) > 0;
    if p1_dom == p2_dom
        % both parents are non-dominated
        if SF(p1) > SF(p2)
            p = p1;
        else
            p = p2;
        end
    elseif ~p1_dom
        p = p1;
    else
        p = p2;
    end
    if rand < p_c
        % find second parent within mating threshold
        searchperm = randperm(N);
        found = false;
        for j = searchperm
            if norm(F(:,i) - F(:,j)) < sigma_mate
                pp = j;
                found = true;
                break;
            end
        end
        if ~found
            % no parent found within mating threshold
            pp = randi(N);
        end
        % intermediate crossover
        alpha = rand;
        X_new(:, i) = alpha * X(:, p) + (1 - alpha) * X(:, pp);
    else
        X_new(:, i) = X(:, i);
    end
end
end

function X_new = mutation(X, p_m, lim, beta)
% Perform mutation
% X: population to mutate
% p_m: mutation probability
% lim: limits of the parameters
dim = size(X, 1);
N = size(X, 2);
X_new = X;
for i = 1:N
    for j = 1:dim
        if rand < p_m
            r = rand;
            if r > 0.5
                delta_x = (lim(j,2) - X(j,i)) * (2 * (r - 0.5))^beta;
            else
                delta_x = (lim(j,1) - X(j,i)) * (2 * (0.5 - r))^beta;
            end
            X_new(j,i) = X(j,i) + delta_x;
        end
    end
end
end
