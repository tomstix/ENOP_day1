clear; clc; close all;

max_iterations = 250; % maximum number of iterations for the optimization

chi0 = 0.7298; % central value for chi
numchi = 5; % number of values to test for chi
delta_chi = 0.5; % spacing between values of chi
chis = linspace(chi0 - delta_chi/2, chi0 + delta_chi/2, numchi);

c1_0 = 2.05; % central value for c1
numcs = 5; % number of values to test for c1 and c2
delta_c = 3; % spacing between values of c1 and c2
c1s = linspace(c1_0 - delta_c/2, c1_0 + delta_c/2, numcs);
c2s = 4.1-c1s;

ns = 50:100:250; % number of particles

num_runs = 100; % number of runs for each experiment. Higher values will give smoother curves

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

    % do the experiments for the different values of chi
    figure
    hold on
    for c = 1:width(chis)
        chi = chis(c);
        c1 = c1s(1);
        c2 = c2s(1);
        num_particles = ns(1);
        f_hists = zeros(num_runs, max_iterations);
        for run = 1:num_runs
            [f_best, g, g_hist, f_hist, k] = particle_swarm_optimization(fn, dim, num_particles, lim, max_iterations, chi, c1, c2, false);
            f_hists(run, :) = f_hist;
        end
        f_hist = mean(f_hists);
        plot(1:k, f_hist(1:k), 'DisplayName', ['\chi = ', num2str(chi), ', Final value: ', num2str(f_hist(k))])
        xscale log
        % yscale log
        xlim([1, k])
    end
    grid on
    title(sprintf('%s (dim = %d, c1 = %0.2f, c2 = %0.2f, particles: %d)', name, dim, c1, c2, num_particles))
    legend
    hold off
    if savefigs
        saveas(gcf, sprintf('figures/%s_chi.png', name))
    end
    
    % do the experiments for the different values of c1 and c2
    figure
    hold on
    for c = 1:width(c1s)
        chi = chis(1);
        c1 = c1s(c);
        c2 = c2s(c);
        num_particles = ns(1);
        f_hists = zeros(num_runs, max_iterations);
        for run = 1:num_runs
            [f_best, g, g_hist, f_hist, k] = particle_swarm_optimization(fn, dim, num_particles, lim, max_iterations, chi, c1, c2, false);
            f_hists(run, :) = f_hist;
        end
        f_hist = mean(f_hists);
        plot(1:k, f_hist(1:k), 'DisplayName', sprintf('c1 = %0.2f, c2 = %0.2f, Final value: %0.4e', c1, c2, f_hist(k)))
        xscale log
        % yscale log
        xlim([1, k])
    end
    grid on
    title(sprintf('%s (dim = %d, \\chi = %s, particles: %s)', name, dim, num2str(chi), num2str(num_particles)))
    legend
    hold off
    if savefigs
        saveas(gcf, sprintf('figures/%s_c.png', name))
    end

    % do the experiments for the different number of particles
    figure
    hold on
    for c = 1:width(ns)
        chi = chis(1);
        c1 = c1s(1);
        c2 = c2s(1);
        num_particles = ns(c);
        f_hists = zeros(num_runs, max_iterations);
        for run = 1:num_runs
            [f_best, g, g_hist, f_hist, k] = particle_swarm_optimization(fn, dim, num_particles, lim, max_iterations, chi, c1, c2, false);
            f_hists(run, :) = f_hist;
        end
        f_hist = mean(f_hists);
        plot(1:k, f_hist(1:k), 'DisplayName', sprintf('particles = %d, Final value: %0.4e', num_particles, f_hist(k)))
        xscale log
        % yscale log
        xlim([1, k])
    end
    grid on
    title(sprintf('%s (dim = %d, \\chi = %s, c1 = %0.2f, c2 = %0.2f)', name, dim, num2str(chi), c1, c2))
    legend
    hold off
    if savefigs
        saveas(gcf, sprintf('figures/%s_n.png', name))
    end
end