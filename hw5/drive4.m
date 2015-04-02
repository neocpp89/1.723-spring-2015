clear all;
close all;
clc;

Ns = [10, 20, 50, 100, 200];
vars = [0.1, 2, 5];
for j=1:numel(vars)
h = figure;
for i=1:numel(Ns);
    [ Cm, ta, p ] = integrate_scalar_tracer_equation_2d( Ns(i), Ns(i), vars(j), 1, 1, 1 );
    Lx = 1;
    hx = Lx / Ns(i);
    x = linspace(0+hx/2.0, Lx-hx/2.0, Ns(i));
    plot(x, diag(p) - p(end/2,end/2)); hold on;
    ll{i} = sprintf('N = %d', Ns(i));
end
legend(ll);
set(h, 'units', 'inches', 'position', [1 1 3 3])
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
fname = sprintf('figs/pressure_%d.png', 10*vars(j));
title('Pressure along diagonal');
xlabel('Coordinate');
ylabel('Pressure');
print(fname, '-dpng');

end
