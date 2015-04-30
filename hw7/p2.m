% calculate and plot errors for problem 2
clear all; close all; clc;

powers = 3:12;
errs = zeros(numel(powers), 1);

for i=1:numel(powers)
    np = 2^powers(i);
    h = 2 * pi / np;
    D = p2_fdmatrix(np);
    D = D ./ (h^2);
    x = linspace(h,2*pi,np)';
    ufn = @(xx) exp(sin(xx).^2);
    ddufn = @(xx) exp(sin(xx).^2).*(2 - 4*(sin(xx).^4));
    u = ufn(x);
    ddu = D*u;
    % figure; hold all; 
    % plot(x, ddu);
    % plot(x, ddufn(x));
    % plot(x, ddu - ddufn(x));
    % hold off;
    errs(i) = sqrt(sum((ddufn(x) - ddu).^2)) / np;
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
title('L_2 error of second derivative of e^{sin^2{x}}');
xlabel('N');
ylabel('L_2 error');
ax = gca; 
ax.XTick = 10.^[1:4];
grid on;
print('p2err.pdf', '-dpdf');

clear all; close all; clc;

powers = 3:12;
errs = zeros(numel(powers), 1);

for i=1:numel(powers)
    np = 2^powers(i);
    h = 2 * pi / np;
    D = p2_fdmatrix(np);
    D = D ./ (h^2);
    x = linspace(h,2*pi,np)';
    ufn = @(xx) exp(sin(xx) .* abs(sin(xx)));
    ddufn = @(xx) exp(sin(xx) .* abs(sin(xx))).*(2*cos(2*xx).*sign(sin(xx)) + (sin(2*xx).^2));
    u = ufn(x);
    ddu = D*u;
    % figure; hold all; 
    % plot(x, ddu);
    % plot(x, ddufn(x));
    % plot(x, ddu - ddufn(x));
    % hold off;
    errs(i) = sqrt(sum((ddufn(x) - ddu).^2)) / np;
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
title('L_2 error of second derivative of e^{sin{x} |sin{x}|}');
xlabel('N');
ylabel('L_2 error');
ax = gca; 
ax.XTick = 10.^[1:4];
grid on;
print('p2errabs.pdf', '-dpdf');

