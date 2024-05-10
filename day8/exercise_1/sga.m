function [xbest,fbest] = sga(fun,N,lb,ub,pc,pm,L,k_max)
pm0 = pm;
P = initialize_population(length(lb),N,L);
[F,xbest,fbest] = evaluate(fun,P,lb,ub);
while k < k_max
    Ptemp = selection(P,F);
    P = crossover(Ptemp,pc);
    P = mutation(P,pm);
    P{1} = xbest;
    [F,xbest,fbest] = evaluate(fun,P,lb,ub,fbest);
    pm = adjust_mutation_rate(P,k,k_max,pm,pm0);
    k = k+1;
end
end

function P = initialize_population(n,N,L)
for j = 1:N
    P{j} = floor(2*rand(1,n*L));
end
end

function [F,xbest,fbest] = evaluate(fun,P,lb,ub,varargin)
j0 = 1; fbest = Inf;
if nargin>4
    j0 = 2;
    xbest = P{1};
    fbest = varargin{1};
    F(1) = fbest;
end
for j = j0:N
    F(j) = feval(decode_bin(P{j},lb,ub));
    if F(j) < fbest
        xbest = P{j};
        fbest = F(j);
    end
end
end

function R = selection(P,F)
N = length(P);
for j = 1:2*N
    c = ceil(N*rand(1,2));
    if F(c(1)) < F(c(2))
        R{j} = P{c(1)};
    else
        R{j} = P{c(2)};
    end
end
end

function P = mutation(P,pm)
N = length(P);
n = length(P{1});
for j = 1:N
    for k = 1:n
        if rand() < pm
            P{j}(k) = 1 - P{j}(k);
        end
    end
end
end

function R = crossover(P,pc)
N = length(P);
n = length(P{1});
for j = 1:N/2
    if rand() < pc
        c = floor(2*rand(1,n));
        for k = 1:n
            R{j}(k) = P{2*j-c(k)}(k);
        end
    else
        R{j} = P{2*j};
    end
end
end

function pm = adjust_mutation_rate(P,k,k_max,pm,pm0)
S = population_diversity(P);
if k < k_max/2
    if S < 0.1
        pm = pm*1.3;
    else pm = pm/1.2;
    end
else
    pm = pm0*(2*(k_max-k)/k_max)^2;
end
end

function S = population_diversity(P)
N = length(P);
n = length(P{1});
for j = 1:n
    for k = 1:N
        L(k) = P{k}(j);
    end
    R(j) = std(L);
end
S = sum(R)/n;
end