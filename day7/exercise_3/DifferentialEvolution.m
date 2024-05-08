clear
clc

% Variable that decides which function is used
functionSelect = 3;
% Valid values are 1, 2 and 3, which function is which can be seen below

% Defining the functions
switch functionSelect
    case 1
        funck = @rosenbrock;
        n = 2; % Can be 2
        range = [-2 2];
    case 2
        funck = @fp;
        n = 2; % Can be 1 or 2
        range = [-2 2];
    case 3
        %funck = @(x) -20 * exp(-0.2 * sqrt(1 / size(x,2) * sum(x(1:end).^2))) - exp(1 / size(x,2) * sum(cos(2 * pi .* x(1:end)))) + 21;
        funck = @auckley;
        n = 2; % Can be 1, 2 or 3
        range = [-10 10];
end

% Plot setup
if n == 2
    figure
    contour = fcontour(@(x,y)[funck([x; y])], range, "LevelStep", 2);
    xlim(range)
    ylim(range)
    hold on
    ax = gca;
    drawnow
end

% Algorithm setup
maxIterations = 100;
iterations = 1;
F = 0.8;
CR = 0.9;
N = 10; % Amount of individuals
bestIndividual = zeros(n, 1);
x = zeros(n,N); % Randomly generated individuals
for i = 1:N
    x(:,i) = (range(2) - range(1)).*rand(n,1) + range(1);
end

% Evaluate starting individuals
previousBest = zeros(n, 1);
currentBest = zeros(n, 1);
eps = 1e-6; % This sets the break condition for the optimization
improvement = eps * 2;
improvementHistory = zeros(n, maxIterations);

for i = 1:N
    if i == 1
        currentBest = x(:,1);
        previousBest = x(:,2);
    elseif funck(x(:,i)) < funck(currentBest)
        previousBest = currentBest;
        currentBest = x(:,i);
    end
end

% Algorithm execution
while iterations < maxIterations & (improvement > eps || improvement == 0)
    for i = 1:N
        % Pick agents
        possibleAgents = 1:N;
        possibleAgents(i) = [];
        agentNumbers = randperm(N-1,3);
        agentA = x(:,possibleAgents(agentNumbers(1)));
        agentB = x(:,possibleAgents(agentNumbers(2)));
        agentC = x(:,possibleAgents(agentNumbers(3)));
        agentX = x(:,i);
        
        % Choose p
        p = randi(n);
        
        % Find potential new position for individual x
        y = zeros(n, 1);
        for j = 1:n
            r = rand;
            if j == p || r < CR
                y(j) = agentA(j) + F * (agentB(j) - agentC(j));
            else
                y(j) = agentX(j);
            end
        end
        if funck(y) < funck(agentX)
            agentX = y;
        end
        x(:,i) = agentX;
    end

    % Evaluate individuals
    for i = 1:N
        if i == 1
            bestIndividual = x(:,1);
        elseif funck(x(:,i)) < funck(bestIndividual)
            bestIndividual = x(:,i);
        end
    end
    
    if funck(bestIndividual) < funck(currentBest)
        previousBest = currentBest;
        currentBest = bestIndividual;
        improvement = abs(funck(previousBest) - funck(currentBest));
    else
        improvement = 0;
    end
    improvementHistory(:,iterations) = improvement;
    
    if n == 2
        scatter(x(1,:), x(2,:), 20, 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1);
        if iterations > 1
            delete(pl);
        end
        pl = scatter(x(1,:), x(2,:), 'filled', 'MarkerFaceColor', [1 0 0]);
        drawnow
        % pause(0.05);
    end

    iterations = iterations + 1;

    % Display information
    disp(['Iteration: ', num2str(iterations), '   CurrentBest: ', num2str(currentBest(1)), ', ', num2str(currentBest(2))]);
    disp(['Function Value: ', num2str(funck(currentBest)), '   Improvement: ', num2str(improvement)])
end

% Display final information
hold off
figure
plot(1:1:maxIterations, improvementHistory);
grid on