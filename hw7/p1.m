% calculate and plot errors for problem 1
clear all; close all; clc;

powers = 3:12;
errs = zeros(numel(powers), 1);

for i=1:numel(powers)
    np = 2^powers(i);
    D = p1_fdmatrix(np);
    h = 1 / (np-1);
    D = D ./ h;
    x = linspace(0,1,np)';
    ufn = @(xx) exp(sin(2*pi*xx));
    dufn = @(xx) 2*pi*cos(2*pi*xx).*exp(sin(2*pi*xx));
    u = ufn(x);
    du = D*u;
    % figure; hold all; 
    % plot(x, du);
    % plot(x, dufn(x));
    % plot(x, du - dufn(x));
    % hold off;
    errs(i) = sqrt(sum((dufn(x) - du).^2)) / np;
end

h = figure;
set(h, 'units', 'inches', 'position', [1 1 4 4])
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
loglog(2.^powers, errs);
title('L_2 error of derivative of e^{sin{2 \pi x}}');
xlabel('N');
ylabel('L_2 error');
ax = gca; 
ax.XTick = 10.^[1:4];
grid on;
print('p1err.pdf', '-dpdf');

