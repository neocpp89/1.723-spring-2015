% solve problem 1
clear all; close all; clc;

% common
dt = 4e-5;
tf = 9;
omega = 1;
Pe = 10000;
c0 = @(x) exp(-200 * (x - 0.7).^2);
grids = [128, 256];

% part 1
for i=1:numel(grids);
    N = grids(i);
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
        u = (up - um) * (theta <= pi) + um;
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
    s = toc;
    fprintf('N = %d took %g seconds with FD.\n', N, s);
    legend('show', 'Location', 'EastOutside');
    title(sprintf('Concentration Profile with Standard FD (N = %d)', N));
    xlabel('x');
    ylabel('c');
    xlim([0, 2*pi]);
    ylim([-0.2, 1]);
    print(sprintf('p2a_%d.pdf', N), '-dpdf');
end

% part 2
for i=1:numel(grids);
    N = grids(i);
    dx = 2*pi / N;
    x = linspace(0, 2*pi - dx, N)';
    c = c0(x);

    [ADD, BDD] = cfd24(N, dx);
    [AD, BD] = cfd14(N, dx);

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
        u = (up - um) * (theta <= pi) + um;

        dcdiff = -(1 / Pe) * (ADD \ (BDD * c));
        dcconv = u .* (AD \ (BD * c));
        dc = -dt * (dcconv + dcdiff);

        c = c + dc;
        t = t + dt;

        if (t > times(j))
            plot(x, c, 'DisplayName', sprintf('t = %f', t));
            j = j + 1;
            if (j > numel(times))
                break;
            end
        end
    end
    s = toc;
    fprintf('N = %d took %g seconds with CFD.\n', N, s);
    legend('show', 'Location', 'EastOutside');
    title(sprintf('Concentration Profile with Compact FD (N = %d)', N));
    xlabel('x');
    ylabel('c');
    xlim([0, 2*pi]);
    ylim([-0.2, 1]);
    print(sprintf('p2b_%d.pdf', N), '-dpdf');
end

% part 3
for i=1:numel(grids);
    N = grids(i);
    dx = 2*pi / N;
    x = linspace(0, 2*pi - dx, N)';
    c = c0(x);
    k1 = 1i*[0:N/2-1 0 -N/2+1:-1]';
    k2 = 1i*[0:N/2 -N/2+1:-1]';
    k2 = k2.*k2;

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
        u = (up - um) * (theta <= pi) + um;

        c_tilde = fft(c);

        dcdiff = -(1 / Pe) * real(ifft(k2 .* c_tilde));
        dcconv = u .* real(ifft(k1 .* c_tilde));
        dc = -dt * (dcconv + dcdiff);

        c = c + dc;
        t = t + dt;

        if (t > times(j))
            plot(x, c, 'DisplayName', sprintf('t = %f', t));
            j = j + 1;
            if (j > numel(times))
                break;
            end
        end
    end
    s = toc;
    fprintf('N = %d took %g seconds with FFT.\n', N, s);
    legend('show', 'Location', 'EastOutside');
    title(sprintf('Concentration Profile with FFT (N = %d)', N));
    xlabel('x');
    ylabel('c');
    xlim([0, 2*pi]);
    ylim([-0.2, 1]);
    print(sprintf('p2c_%d.pdf', N), '-dpdf');
end
