% solve problem 1
clear all; close all; clc;

% common
dt = 4e-5;
tf = 9;
omega = 1;
Pe = 10000;
c0 = @(x) exp(-200 * (x - 0.7).^2);

% part 1
N = 256;
dx = 2*pi / N;
x = linspace(0, 2*pi - dx, N)';
c = c0(x);

% DF = fd14_forward_periodic(N, dx);
% DB = fd14_backward_periodic(N, dx);
DDC = fd24_central_periodic(N, dx);
DC = fd14_central_periodic(N, dx);

up = 2 * sin(x / 2);
um = -up;

t = 0;
times = [0 1.73 4.55 8.01];
j = 1;

h = figure;
set(h, 'units', 'inches', 'position', [1 1 6 4]);
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
hold all;
tic;
while (t < tf)
    theta = mod(omega * t, 2*pi);
    % u = (up - um) * (theta <= pi) + um;
    % dcconv = u.* ((DB * c) .* (u >= 0) + (DF * c) .* (u < 0));

    % cppp = [c(4:end); c(1); c(2); c(3);];
    % cpp = [c(3:end); c(1); c(2);];
    % cp = [c(2:end); c(1)];
    % cm = [c(end); c(1:end-1)];
    % cmm = [c(end-1); c(end); c(1:end-2)];
    % cmmm = [c(end-2); c(end-1); c(end); c(1:end-3)];

    % dcconv = (u > 0) .* u .* (c - cm) / (dx) + (u <= 0) .* u .* (cp - c) / (dx);
    % dcconv = u .* ((u > 0) .* (3*c - 4*cm + 1*cmm) / (2*dx) + (u <= 0) .* (-1*cpp + 4*cp - 3*c) / (2*dx));
    % dcconv = u .* ((u > 0) .* (2*cp + 3*c - 6*cm + 1*cmm) / (6*dx) + (u <= 0) .* (-2*cm - 3*c + 6*cp - 1*cpp) / (6*dx));

    Adiff = -(1 / Pe) * DDC;
    dcdiff = Adiff * c;
    dcconv = u .* (DC * c);
    dc = -dt * (dcconv + dcdiff);

    c = c + dc;
    t = t + dt;

    % plot(x, c);
    % title(sprintf('t = %f', t));
    % drawnow;

    if (t > times(j))
        plot(x, c, 'DisplayName', sprintf('t = %f', t));
        j = j + 1;
        if (j > numel(times))
            break;
        end
    end
end
toc;
legend('show', 'Location', 'EastOutside');
title(sprintf('Concentration Profile with Standard FD (N = %d)', N));
xlabel('x');
ylabel('c');
xlim([0, 2*pi]);
ylim([-0.2, 1]);
print('p2a.pdf', '-dpdf');

