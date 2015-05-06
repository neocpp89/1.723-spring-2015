% solve problem 1
clear all; close all; clc;

dx = 0.01;
kmax = pi / dx;
k = linspace(0, kmax, 1000);

g2 = @(k,dx) (sin(k * dx) / dx);
g4 = @(k,dx) ((8*sin(k*dx) - sin(2*k*dx))/(6*dx));
ginf = @(k) k;

h = figure;
set(h, 'units', 'inches', 'position', [1 1 4 4])
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
hold all;
plot(k, g2(k,dx), 'DisplayName', 'g_2');
plot(k, g4(k,dx), 'DisplayName', 'g_4');
plot(k, ginf(k), 'DisplayName', 'g_\infty');
legend('show', 'Location', 'NorthWest');
title('Imaginary parts of differentiation functions');
xlabel('k');
ylabel('Im(g)');
xlim([0, kmax]);
ylim([0, kmax]);
print('p1.pdf', '-dpdf');
