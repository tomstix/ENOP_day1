function [f_best, f_hist, p_best, P, pm, evals] = fga(f, dim, lim, num_individuals, num_generations, pc, pm0, beta, dispfun)
k = 0; % initialize the iteration counter
evals = 0; % initialize the number of function evaluations
pm = pm0; % initialize the mutation probability
f_hist = zeros(1,num_generations); % initialize the history of the best fitness values
P = initialize_population(num_individuals, dim, lim);
F = evaluate_population(P, f);
evals = evals + num_individuals;
while k < num_generations
    k = k+1;
    pairs = selection(P,F);
    P_new = crossover(pairs,pc);
    P_new = mutation(P_new,pm,beta,lim);
    P = P_new;
    F = evaluate_population(P, f);
    f_hist(k) = min(F);
    evals = evals + num_individuals;
    pm = adjust_mutation_rate(P,k,num_generations,pm,pm0);
    if nargin == 9
        dispfun(k, P, F, pm, evals);
    end
end
[f_best, idx] = min(F);
p_best = P(:,idx);
end

function P = initialize_population(n, dim, lim)
P = lim(1) + (lim(2) - lim(1)) * rand(dim, n);
end

function F = evaluate_population(P, f)
F = f(P);
end

function R = selection(P,F)
N = length(P);
dim = height(P);
R = zeros(dim,2*N);
for j = 1:2*N
    c = ceil(N*rand(1,2));
    if F(c(1)) < F(c(2))
        R(:,j) = P(:,c(1));
    else
        R(:,j) = P(:,c(2));
    end
end
end

function R = crossover(P,pc)
N = length(P);
dim = height(P);
R = zeros(dim,N/2);
for j = 1:N/2
    if rand() < pc
        % intermediate crossover
        alpha = rand();
        R(:,j) = alpha*P(:,2*j) + (1-alpha)*P(:,2*j-1);
    else
        % copy one parent
        R(:,j) = P(:,2*j);
    end
end
end

function R = mutation(P,pm,beta,lim)
N = length(P);
dim = height(P);
for j = 1:N
    if rand() < pm
        for k = 1:dim
            r = rand();
            if r > 0.5
                delta_x = (lim(2) - P(k,j))*(2*(r - 0.5))^beta;
            else
                delta_x = (lim(1) - P(k,j))*(2*(0.5-r))^beta;
            end
            P(k,j) = P(k,j) + delta_x;
        end
    end
end
R = P;
end

function pm = adjust_mutation_rate(P,k,k_max,pm,pm0)
S = population_diversity(P);
if k < k_max/2
    if S < 0.1
        pm = pm*1.3;
    else
        pm = pm/1.2;
    end
else
    pm = pm0*(2*(k_max-k)/k_max)^2;
end
end

function S = population_diversity(P)
% calculate the standard deviation of the distances from the centroid to the points
N = length(P);
cent = mean(P,2);
dist = P - cent;
norms = zeros(1,N);
for j = 1:N
    norms(j) = norm(dist(:,j));
end
S = std(norms);
end