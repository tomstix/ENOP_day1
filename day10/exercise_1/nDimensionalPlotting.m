function nDimensionalPlotting(sampleSet)
    n = size(sampleSet, 2); % Dimension of sampleSet
    subSpaceCount = (n*(n-1)/2); % Amount of subspaces
    subSpaces = nchoosek(1:n,2); % A matrix containing all the subspaces
    
    plotMatrixSize = factor(subSpaceCount);
    if size(plotMatrixSize, 2) == 1
        plotMatrixSize = [1, plotMatrixSize(end)];
    elseif size(plotMatrixSize, 2) > 2
        plotMatrixSize = [prod(plotMatrixSize(1:end-1)), plotMatrixSize(end)];
    end

    figure
    for i = 1:subSpaceCount
        subplot(plotMatrixSize(1), plotMatrixSize(2), i)
        scatter(sampleSet(:,subSpaces(i,1)), sampleSet(:,subSpaces(i,2)))
        xlabel(['x_{', num2str(subSpaces(i,1)), '}'], 'FontWeight', 'bold')
        ylabel(['x_{', num2str(subSpaces(i,2)), '}'], 'FontWeight', 'bold')
    end
end