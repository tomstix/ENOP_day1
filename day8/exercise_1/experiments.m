clear; clc; close all;

num_individuals = 10:30:100; % number of individuals
num_generations = 200; % number of generations
pc0 = 0.8; % crossover probability
pm0_0 = 0.1; % initial mutation probability
beta = 3; % mutation parameter

num_pc = 5; % number of different values of pc to try
num_pm0 = 5; % number of different values of pm0 to try
delta_pc = 0.2; % spacing between values of pc
delta_pm0 = 0.2; % spacing between values of pm0
pcs = linspace(pc0 - delta_pc/2, pc0 + delta_pc/2, num_pc);
pm0s = linspace(pm0_0 - delta_pm0/2, pm0_0 + delta_pm0/2, num_pm0);

num_runs = 50; % number of runs for each experiment. Higher values will give smoother curves

savefigs = false; % set to true to save the figures

for function_select = 1:3 % loop over the functions
    switch function_select
        case 1
            fn = @rosenbrock;
            lim = [-2 2];
            dims = 2;
            name = 'Rosenbrock';
        case 2
            fn = @fp;
            lim = [-2 2];
            dims = [1,2];
            name = 'f_p';
        case 3
            fn = @auckley;
            lim = [-10 10];
            dims = [1,2,3];
            name = 'Auckley';
    end
    dim = 2; % only do the 2D case to keep the number of plots manageable
    
    % do the experiments for the different values of num_individuals
    figure
    hold on
    for c = 1:width(num_individuals)
        individuals = num_individuals(c);
        pc = pc0;
        pm0 = pm0_0;
        f_hists = zeros(num_runs, num_generations);
        for run = 1:num_runs
            [f_best, f_hist, p_best, P, pm, evals] = fga(fn, dim, lim, individuals, num_generations, pc, pm0, beta);
            f_hists(run, :) = f_hist;
        end
        f_hist = mean(f_hists);
        plot(1:num_generations, f_hist(1:num_generations), 'DisplayName', ['Individuals = ', num2str(individuals), ', Final value: ', num2str(f_best)])
        xscale log
        % yscale log
        xlim([1, num_generations])
    end
    grid on
    title(sprintf('%s (dim = %d, pc = %0.2f, pm0 = %0.2f, beta = %0.2f)', name, dim, pc, pm0, beta))
    legend
    hold off
    if savefigs
        saveas(gcf, sprintf('figures/%s_ind.png', name))
    end
    
    % do the experiments for the different values of pc
    figure
    hold on
    for c = 1:width(pcs)
        pc = pcs(c);
        individuals = num_individuals(1);
        pm0 = pm0_0;
        f_hists = zeros(num_runs, num_generations);
        for run = 1:num_runs
            [f_best, f_hist, p_best, P, pm, evals] = fga(fn, dim, lim, individuals, num_generations, pc, pm0, beta);
            f_hists(run, :) = f_hist;
        end
        f_hist = mean(f_hists);
        plot(1:num_generations, f_hist(1:num_generations), 'DisplayName', ['pc = ', num2str(pc), ', Final value: ', num2str(f_best)])
        xscale log
        % yscale log
        xlim([1, num_generations])
    end
    grid on
    title(sprintf('%s (dim = %d, individuals = %d, pm0 = %0.2f, beta = %0.2f)', name, dim, individuals, pm0, beta))
    legend
    hold off
    if savefigs
        saveas(gcf, sprintf('figures/%s_pc.png', name))
    end

    % do the experiments for the different values of pm0
    figure
    hold on
    for c = 1:width(pm0s)
        pm0 = pm0s(c);
        individuals = num_individuals(1);
        pc = pc0;
        f_hists = zeros(num_runs, num_generations);
        for run = 1:num_runs
            [f_best, f_hist, p_best, P, pm, evals] = fga(fn, dim, lim, individuals, num_generations, pc, pm0, beta);
            f_hists(run, :) = f_hist;
        end
        f_hist = mean(f_hists);
        plot(1:num_generations, f_hist(1:num_generations), 'DisplayName', ['pm0 = ', num2str(pm0), ', Final value: ', num2str(f_best)])
        xscale log
        % yscale log
        xlim([1, num_generations])
    end
    grid on
    title(sprintf('%s (dim = %d, individuals = %d, pc = %0.2f, beta = %0.2f)', name, dim, individuals, pc, beta))
    legend
    hold off
    if savefigs
        saveas(gcf, sprintf('figures/%s_pm0.png', name))
    end
end