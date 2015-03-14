clear all;
close all;
clc;

x = linspace(0,1,200);
times = [1e-3, 1e-1, 1, 2]
legend_labels = arrayfun(@(t) sprintf('Time: %g', t), times, 'uniform', false);
Pe_range = [1, 10, 100]

for Pe=Pe_range
    h = figure;
    hold off;
    for t=times
        plot(x, 0.5*erfc((x-t)/(2*sqrt(t/Pe))), 'LineWidth', 2);
        hold all;
    end
    axis equal tight;
    set(h, 'units', 'inches', 'position', [1 1 3 3])
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    title(sprintf('Linear Tracer Transport - Pe = %g', Pe), 'interpreter', 'latex')
    xlabel('$x_D$', 'interpreter', 'latex')
    ylabel('$c_D$', 'interpreter', 'latex')
    xlim([0 1])
    ylim([0 1])
    legend(legend_labels)
    fname = sprintf('figs/analytic_%g.pdf', Pe)
    print(fname, '-dpdf')
end
