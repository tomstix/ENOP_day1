function sampleSet = HypercubeSampling(N, n)
    sampleSet = zeros(N,n);
    binSize = 1/N; % Size of the bins in one dimension

    for i = 1:n
        order = randperm(N);
        for j = 1:N
            binUpperBound = order(j) / N;
            binLowerBound = binUpperBound - binSize;
            sampleSet(j,i) = (binUpperBound-binLowerBound) * rand + binLowerBound;
        end
    end
end