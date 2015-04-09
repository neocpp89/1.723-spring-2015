clear all;
close all;
clc;

alpha = 1;
Nx = 128;
Nz = 128;
t_wanted = [1, 4, 10, 15];
tf = 20;
drawdelay = 500;
dt = 0.001;

for Ra=[1000, 2000, 4000]
    [cm, tm, ct, tt] = integrate_density(alpha, Nx, Nz, Ra, dt, tf, drawdelay, t_wanted);
    Lz = 1;
    hz = Lz / Nz;
    z = linspace(0+hz/2.0, Lz-hz/2.0, Nz);
    h = figure;
    set(h, 'units', 'inches', 'position', [1 1 5 5])
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    hold on;
    for j=1:size(cm,3)
        cmz = mean(cm(:,:,j), 2);
        plot(z, cmz);
    end
    title('Horizontally Averaged Concentration');
    xlabel('Z');
    ylabel('C');
    fname = sprintf('figs/Ra_%d_cbar.png', Ra);
    legend(t_wanted, 'location', 'northeast');
    print(fname, '-dpng');
    hold off;

    h = figure;
    set(h, 'units', 'inches', 'position', [1 1 5 5])
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    tt = [0, tt];
    ct = [0, ct];
    plot(tt, ct);
    title('Total Concentration');
    xlabel('T');
    ylabel('C');
    fname = sprintf('figs/Ra_%d_ctot.png', Ra);
    print(fname, '-dpng');

end
