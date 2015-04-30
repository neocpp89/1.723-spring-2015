clear all;
close all;
clc;

N = 40;
Pe = [100, 200, 500, 1000, 2000];
vars = [0.1, 2, 5];
for j=1:numel(vars)
h = figure;
for i=1:numel(Pe);
    [ Cm, ta, p ] = integrate_scalar_tracer_equation_2d(N, N, vars(j), Pe(i), 0.0001, 2);
    plot(ta, Cm); hold on;
    ll{i} = sprintf('Pe = %d', Pe(i));
end
legend(ll, 'Location', 'best');
set(h, 'units', 'inches', 'position', [1 1 3 3])
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition', [0 0 pos(3) pos(4)]);
fname = sprintf('figs/breakthrough_%d.png', 10*vars(j));
title('Breakthough Curve');
xlabel('Time');
ylabel('Concentration');
xlim([0, 2]);
ylim([0, 1]);
print(fname, '-dpng');

end
