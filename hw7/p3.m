% solve problem 3
clear all; close all; clc;

% want polynomial p(x) = ax^3 + bx^2 + cx + d

% part 2
N = 4;
x = linspace(-1, 1, N);
y = sin(x);
[p,r,s] = lagrangepoly(x, y);
xx = -1:0.1:1;
yy = polyval(p, xx);

h = figure;
set(h, 'units', 'inches', 'position', [1 1 4 4])
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
hold all;
plot(xx, sin(xx), '-ob', 'DisplayName', 'Sampled');
plot(xx, yy, '-xr', 'DisplayName', 'Intepolated');
legend('show', 'Location', 'SouthEast');
title('Sampled vs Interpolated (N=4) sin(x)');
xlabel('x');
ylabel('y');
print('p3n4.pdf', '-dpdf');

% part 3
fn = @(g) 1 ./ (1 + 16*(g.^2));
h = figure;
for N=[4 20]
    x = linspace(-1, 1, N);
    y = fn(x);
    [p,r,s] = lagrangepoly(x, y);
    xx = -1:0.1:1;
    yy = polyval(p, xx);
    set(h, 'units', 'inches', 'position', [1 1 4 4])
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    hold all;
    if (N == 4)
        plot(xx, fn(xx), '-o', 'DisplayName', 'Sampled');
    end
    plot(xx, yy, '-x', 'DisplayName', sprintf('Interpolated (N=%d)', N));
end
legend('show', 'Location', 'NorthEast');
title('Sampled vs Interpolated 1/(1+16x^2)');
xlabel('x');
ylabel('y');
print('p3n20.pdf', '-dpdf');

% part 4
fn = @(g) 1 ./ (1 + 16*(g.^2));
h = figure;
for N=[4 20]
    n = 1:N;
    x = cos(pi*(n-1)/(N - 1));
    y = fn(x);
    [p,r,s] = lagrangepoly(x, y);
    xx = -1:0.1:1;
    yy = polyval(p, xx);
    set(h, 'units', 'inches', 'position', [1 1 4 4])
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    hold all;
    if (N == 4)
        plot(xx, fn(xx), '-o', 'DisplayName', 'Sampled');
    end
    plot(xx, yy, '-x', 'DisplayName', sprintf('Nonuniform Interpolated (N=%d)', N));
end
legend('show', 'Location', 'NorthEast');
title('Sampled vs Interpolated 1/(1+16x^2)');
xlabel('x');
ylabel('y');
print('p3n20c.pdf', '-dpdf');
