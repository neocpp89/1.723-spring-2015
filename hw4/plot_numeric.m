clear all;
close all;
clc;

wanted_times = [1e-3, 1e-1, 1, 2];
Pe_range = [1, 10, 100];
N = 100;

analytical_solution = @(x,t,Pe) 0.5*erfc((x-t)/(2*sqrt(t/Pe)));

for dt=[1e-2, 1e-3]
    for Pe=Pe_range
        for theta=[0, 1]
            h = figure;
            hold off;
            [X,C,times] = integrate_scalar_tracer_equation(N, Pe, theta, dt, wanted_times);
            for i=1:numel(times)
                plot(X(:, i), analytical_solution(X(:, i), times(i), Pe), '-', 'LineWidth', 2);
                hold all;
            end
            hold on;
            plot(X, C, 'o');
            axis tight;
            set(h, 'units', 'inches', 'position', [1 1 3.5 3.5]);
            set(h, 'PaperUnits','centimeters');
            set(h, 'Units','centimeters');
            pos=get(h,'Position');
            set(h, 'PaperSize', [pos(3) pos(4)]);
            set(h, 'PaperPositionMode', 'manual');
            set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
            title(sprintf('Linear Tracer Transport\n(Pe = %g, theta = %d, dt = %.1e)', Pe, theta, dt), 'interpreter', 'latex');
            xlabel('$x_D$', 'interpreter', 'latex');
            ylabel('$c_D$', 'interpreter', 'latex');
            % xlim([0 1])
            % ylim([0 max(max(C(:)), 1)])
            xlim([0 1]);
            ylim([0 1]);
            legend_labels = arrayfun(@(t) sprintf('Time: %g', t), times, 'uniform', false);
            legend_labels_n = cellfun(@(t) sprintf('%s - Numeric', t), legend_labels, 'uniform', false);
            legend([legend_labels; legend_labels_n]);
            legend boxoff;
            fname = sprintf('figs/numeric_%d_%d_%d.pdf', Pe, 1000*dt, theta);
            print(fname, '-dpdf');
        end
    end
end
