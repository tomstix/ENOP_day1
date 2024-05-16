clear; clc; close all;

fn_select = 2; % select the function to optimize
switch fn_select
    case 1
        fn = @ft;
        lim = [0.5, 1.0; -2.0, 2.0; -2.0, 2.0];
        dim = 3;
        m = 2;
        plot_limits = [0.4, 1; 1, 5];
        ann_pos = [0.8, 0.8, 0.1, 0.1];
    case 2
        fn = @ff;
        lim = [-2 * ones(8,1), 2 * ones(8,1)];
        dim = 8;
        m = 2;
        plot_limits = [0, 1; 0, 1];
        ann_pos = [0.15, 0.4, 0.1, 0.1];
end

% Parameters
N = 100; % Population size
k_max = 200; % Maximum number of iterations
p_c = 0.8; % Crossover probability
p_m = 0.2; % Mutation probability
n_elite = 40; % Number of elite individuals
alpha_s = 1; % Sharing function shape parameter
update_sigma_mate = @(sigma_share) 3*sigma_share; % Function to update the mating threshold
rand_beta = 1; % Beta for random mutation
break_history = zeros(1, k_max);
break_threshold = 0.02; % Threshold for breaking the loop
pause_time = 0; % Pause time between iterations for plotting

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
% make an annotation with N, n_elite, p_c, p_m, alpha_s, rand_beta
annotation('textbox', ann_pos, 'String', {['N: ', num2str(N)], ['n_{elite}: ', num2str(n_elite)], ['p_c: ', num2str(p_c)], ['p_m: ', num2str(p_m)],  ['\sigma_{mate}: 3*\sigma_{share}'], ['\alpha_s: ', num2str(alpha_s)], ['\beta: ', num2str(rand_beta)]}, 'FitBoxToText', 'on');

nondom_idx = D == 0;
nondom = X(:,nondom_idx);

% Main loop
for k = 1:k_max
    % find the two points with the biggest distance in F
    combs = nchoosek(1:N, 2); % find all combinations of 2 points
    max_dist = 0;
    best_comb = [0, 0];
    for i = 1:size(combs, 1) % iterate through all combinations
        d = norm(F(:, combs(i,1)) - F(:, combs(i,2)));
        if d > max_dist % if the distance is bigger than the current max distance, update the max distance and the best combination
            max_dist = d;
            best_comb = combs(i,:);
        end
    end
    % approximate the length of the pareto front
    F_x = F(:, best_comb(1));
    F_y = F(:, best_comb(2));
    d1 = abs(F_x(2) - F_y(2));
    d2 = abs(F_x(1) - F_y(1));
    d_min = sqrt(d1^2 + d2^2);
    d_max = d1 + d2;
    d_k = (d_min + d_max) / 2;
    % update the sharing function parameters
    sigma_share = N^(1/1-m) * d_k;
    sigma_mate = update_sigma_mate(sigma_share);
    
    % Selection
    X_new = crossover(X, F, D, SF, p_c, sigma_mate);
    X_new = mutation(X_new, p_m, lim, rand_beta);
    X = [X_new, X_elite]; % Add elite individuals to the population
    F = fn(X);  % Evaluate the new population
    D = zeros(1, size(X, 2));
    for i = 1:size(X, 2)
        D(i) = count_dominated_by(F(:,i), F); % count the dominance numbers
    end
    SF = compute_shared_fitnesses(X, D, sigma_share, alpha_s); % Compute the shared fitnesses
    % sort the population by shared fitness
    [~, idx] = sort(SF, "descend");
    X = X(:,idx(1:N));
    F = F(:,idx(1:N));
    D = D(idx(1:N));
    SF = SF(idx(1:N));
    % select the elite individuals using recurrence mode
    X_elite = X;
    while size(X_elite, 2) > n_elite
        F_e = fn(X_elite);
        D_e = zeros(1, size(X_elite, 2));
        for i = 1:size(X_elite, 2)
            D_e(i) = count_dominated_by(F_e(:,i), F_e);
        end
        SF_e = compute_shared_fitnesses(X_elite, D_e, sigma_share, alpha_s);
        [~, min_idx] = min(SF_e);
        X_elite(:,min_idx) = [];
    end
    % Plot
    refreshdata
    title(['Iteration ', num2str(k), ' \sigma_{share} = ', num2str(sigma_share), ' \sigma_{mate} = ', num2str(sigma_mate)]);
    drawnow limitrate;
    
    fprintf("Iteration %d \t d_k: %f \t sigma_share: %f \t sigma_mate: %f\n", k, d_k, sigma_share, sigma_mate)

    % get the new non-dominated solutions
    nondom_idx_new = D == 0;
    nondom_new = X(:,nondom_idx_new);
    sum_new_d = 0;
    F_new = fn(nondom_new);
    F_old = fn(nondom);
    % count how many new non-dominated solutions dominate the old ones
    for i = 1:size(F_new, 2)
        for j = 1:size(F_old, 2)
            if all(F_new(:,i) <= F_old(:,j)) && any(F_new(:,i) < F_old(:,j))
                sum_new_d = sum_new_d + 1;
                break;
            end
        end
    end
    p_k = sum_new_d / size(F_new, 2); % calculate the ratio

    nondom = nondom_new;

    % keep track of the break history to have a moving average
    break_history(k) = p_k;
    if k > 10
        p_k_avg = mean(break_history(k-10:k));
        if p_k_avg < break_threshold
            break; % break if there are almost no new non-dominated solutions
        end
    end
    
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

function X_new = crossover(X, F, D, SF, p_c, sigma_mate)
dim = size(X, 1);
N = size(X, 2);
X_new = zeros(dim, N);
for i = 1:N
    % tournament selection for parent 1
    % select 2 individuals randomly
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
            d = norm(F(:,i) - F(:,j));
            if d ~= 0 && d < sigma_mate
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
            % for each parameter, generate a random delta with non-uniform distribution
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
