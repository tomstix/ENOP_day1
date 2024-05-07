clear
clc

% Variable that decides which function Newton's method is used on
functionSelect = 1;
% Change the above variable and rerun the program to see Newton's method
% used on the 4 different functions.
% Valid values are 1, 2, 3 and 4, which function is which can be seen below

% Defining the functions
switch functionSelect
    case 1
        funck = @(x,y) 2*x.^2 + 3*y.^2 - 3*x.*y + x;
        x_0 = [5 8];
        range = [-2 8];
    case 2
        funck = @(x,y) (1 - x).^2 + 5*(x - y.^2).^2;
        x_0 = [0 0];
        range = [-0.5 1.5];
    case 3
        funck = @(x,y) (x + 2*y).*(1 - 0.9*exp(-0.3*(x - 2.5).^2 - 2*(y - 3.5).^2)).*(1 - 0.9*exp(-(x - 3).^2 - (y - 3).^2));
        x_0 = [4 2];
        range = [1 5];
    case 4
        funck = @(x,y) exp(x/5) + exp(y/3);
        x_0 = [5 8];
        range = [-10 10];
end

% Plot setup
figure
contour = fcontour(funck, range, "LevelStep", 2);
xlim(range)
ylim(range)
hold on
ax = gca;
drawnow

% Patteren search setup
maxSteps = 100;
dim = size(x_0,2);
stepHistory = zeros(maxSteps, dim);
stepHistory(1,:) = x_0;
stepHistory(2,:) = x_0;
commandHistory = strings(maxSteps, 1);
gridSize = 1;
gridReduction = 10;
stepsTaken = 1;
evaluationsMade = 0;
eps = 1e-6; % This sets the break condition for the optimization

% Patteren search execution
while gridSize > eps & stepsTaken < maxSteps
    % Check both neighbours in each dimension
    validNeighbourFound = 0;
    direction = 0;
    for i = 1:dim
        m = zeros(1, dim);
        m(i) = 1;
        fwd = stepHistory(stepsTaken,:) + m*gridSize;
        bwd = stepHistory(stepsTaken,:) - m*gridSize;
        
        % Evaluate both neighbours in this dimension
        % Plot evaluation
        pl = plot(ax, [stepHistory(stepsTaken, 1), fwd(1)], [stepHistory(stepsTaken, 2), fwd(2)], 'g', 'LineWidth',1);
        pause(0.05);
        delete(pl);
        
        evaluationsMade = evaluationsMade + 1;
        if funck(fwd(1), fwd(2)) < funck(stepHistory(stepsTaken + 1,1), stepHistory(stepsTaken + 1,2))
            stepHistory(stepsTaken + 1,:) = fwd;
            stepHistory(stepsTaken + 2,:) = fwd;
            validNeighbourFound = 1;
            direction = i;
            continue
        end
        
        % Plot evaluation
        pl = plot(ax, [stepHistory(stepsTaken, 1), bwd(1)], [stepHistory(stepsTaken, 2), bwd(2)], 'g', 'LineWidth',1);
        pause(0.05);
        delete(pl);

        evaluationsMade = evaluationsMade + 1;
        if funck(bwd(1), bwd(2)) < funck(stepHistory(stepsTaken + 1,1), stepHistory(stepsTaken + 1,2))
            stepHistory(stepsTaken + 1,:) = bwd;
            stepHistory(stepsTaken + 2,:) = bwd;
            validNeighbourFound = 1;
            direction = i + dim;
            continue
        end
        % Plot evaluation
        pl = plot(ax, [stepHistory(stepsTaken, 1), bwd(1)], [stepHistory(stepsTaken, 2), bwd(2)], 'g', 'LineWidth',1);
        pause(0.05);
        delete(pl);
    end
    % If no valid neighbour was found, reduce the grid size
    if validNeighbourFound == 0
        stepHistory(stepsTaken + 1,:) = stepHistory(stepsTaken,:);
        gridSize = gridSize / gridReduction;
    % If a valid neighbour was found, save the step and continue in the same
    % direction
    else
        stepsTaken = stepsTaken + 1;
        visualizePath(ax, stepHistory(1:stepsTaken,:), false);
        
        % Calcuelate the next step location
        dimDirection = 1;
        if direction / dim > 1
            direction = direction - dim;
            dimDirection = -1;
        end
        m = zeros(1, dim);
        m(direction) = dimDirection;

        % Continue in the same direction as long as it gets closer to the
        % minimum
        while stepsTaken < maxSteps
            nextStep = stepHistory(stepsTaken,:) + m*gridSize;

            % Plot evaluation
            pl = plot(ax, [stepHistory(stepsTaken, 1), nextStep(1)], [stepHistory(stepsTaken, 2), nextStep(2)], 'g', 'LineWidth',1);
            pause(0.05);
            delete(pl);
            
            % If the step isn't closer to the minimum, go back to neighbour
            % search
            evaluationsMade = evaluationsMade + 1;
            if funck(nextStep(1), nextStep(2)) > funck(stepHistory(stepsTaken,1), stepHistory(stepsTaken,2))
                break;
            % If the step is closer to the minimum, save the step and
            % continue in the same direction
            else
                stepHistory(stepsTaken + 1,:) = nextStep;
                stepHistory(stepsTaken + 2,:) = nextStep;
                
                stepsTaken = stepsTaken + 1;
                visualizePath(ax, stepHistory(1:stepsTaken,:), false);
            end
        end
    end
end 

% Plot the final path
visualizePath(ax, stepHistory(1:stepsTaken,:), true);