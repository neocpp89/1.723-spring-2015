function [ cm, ta, p ] = integrate( alpha, Nx, Nz, Ra, dt, tf, drawdelay )
%INTEGRATE_SCALAR_TRACER_EQUATION_2D Solves 2D linear scalar transport equation
%using finite volumes.
%  Nx - number of grid cells in x
%  Ny - number of grid cells in x
%  Pe - Peclet number
%  dt - time step size
%  tf - max time
% 
%  Returns Nx x Ny x (# times) matrix xm containing cell point locations,
%          1 x (# times) matrix cm containing top right cell concentrations at different times,
%          1 x (# times) matrix ta containing actual times data is saved at.
%          Nx x Ny matrix containing pressure field.
if (nargin == 6)
    drawdelay = 1000;
end

epsrnd = 0.005;

Lx = alpha;
Lz = 1;
hx = Lx / Nx;
hz = Lz / Nz;
x = linspace(0+hx/2.0, Lx-hx/2.0, Nx);
z = linspace(0+hz/2.0, Lz-hz/2.0, Nz);
[X,Z] = meshgrid(x,z);

lambda = ones(Nz, Nx);
lambda_bar_x = 0.5*(lambda(:, 2:end) + lambda(:, 1:end-1));
lambda_bar_z = 0.5*(lambda(2:end, :) + lambda(1:end-1, :));
lambda_bar_x = (hz / hx) .* lambda_bar_x;
lambda_bar_z = (hx / hz) .* lambda_bar_z;

% Set no flow BCs
lambda_bar_x = [zeros(Nz, 1), lambda_bar_x, zeros(Nz, 1)];
lambda_bar_z = [zeros(1, Nx); lambda_bar_z; zeros(1, Nx)];

b = zeros(Nx,Nz);
bv = reshape(b, [], 1);

% make diagonals
T_x_left = reshape(lambda_bar_x(:, 1:end-1), [], 1);
T_x_right = reshape(lambda_bar_x(:, 2:end), [], 1);
T_z_bottom = reshape(lambda_bar_z(1:end-1, :), [], 1);
T_z_top = reshape(lambda_bar_z(2:end, :), [], 1);

sumT = T_x_right + T_x_left + T_z_bottom + T_z_top;
diags = [-T_x_left, -T_z_bottom, sumT, -T_z_top, -T_x_right];

% assemble matrix
nr = Nx*Nz;
A = spdiags(diags, [Nx;1;0;-1;-Nx], nr, nr)';

assert(A(1,1) == (hx/hz + hz/hx));
assert(A(1,2) == -(hx/hz));
assert(A(1,1+Nx) == -(hz/hx));

% decompose matrix
[L,D,perm] = ldl(A,'vector');

t = 0;
i = 1;
c = zeros(Nz, Nx);
% c(end/3:2*end/3, end/3:2*end/3) = 1;
% c(end/2, end/2) = 1;
% c(1, :) = 1 + 2*0.5*(rand(1,Nx)-0.5);
while (t < tf)
    t = t + dt;
    ta(i) = t;

    c_z_avg = hx*0.5*(c(1:end-1, :) + c(2:end, :));
    c_z = [zeros(1, Nx); c_z_avg; zeros(1, Nx)];
    b = (c_z(1:end-1,:) - c_z(2:end,:)); 

    bv = reshape(b, [], 1);

    % solve for pressure
    clear pv;
    pv(perm,:) = L'\(D\(L\(bv(perm,:))));

    p = reshape(pv, Nz, Nx);

    % show plots and save pressure
    % figure; surf(X,Y,log(perm),'edgecolor','none'); view([0, 0, 1]);
    % figure; surf(X,Y,p,'edgecolor','none'); view([0, 0, 1]);
    % h = figure; plot(x, diag(p) - p(end/2,end/2));
    % set(h, 'units', 'inches', 'position', [1 1 3 3])
    % set(h, 'PaperUnits','centimeters');
    % set(h, 'Units','centimeters');
    % pos=get(h,'Position');
    % set(h, 'PaperSize', [pos(3) pos(4)]);
    % set(h, 'PaperPositionMode', 'manual');
    % set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    % fname = sprintf('figs/pressure_%d_%d.png', 10*var_lnk, Nx);
    % title(sprintf('Pressure along diagonal\n(N = %d)', Nx));
    % xlabel('Coordinate');
    % ylabel('Pressure');
    % print(fname, '-dpng');

    % calculate integrated velocities
    ux_int = lambda_bar_x(:, 2:end-1).*(p(:, 1:end-1, :) - p(:, 2:end));
    uz_int = lambda_bar_z(2:end-1, :).*(p(1:end-1, :) - p(2:end, :));
    uz_int = uz_int + c_z_avg;
    ux = [zeros(Nz, 1), ux_int, zeros(Nz, 1)];
    uz = [zeros(1, Nx); uz_int; zeros(1, Nx)];

    c_top = 1 + 2*epsrnd*(rand(1,Nx)-0.5);
    % c_top = zeros(1, Nx);

    Fx_adv = ux_int.*(ux_int > 0).*c(:, 1:end-1) + ux_int.*(ux_int < 0).*c(:, 2:end);
    Fz_adv = uz_int.*(uz_int > 0).*c(1:end-1, :) + uz_int.*(uz_int < 0).*c(2:end, :);

    Fz_diff = (1 / Ra) * (hx / hz) * (c(1:end-1, :) - c(2:end, :));
    Fx_diff = (1 / Ra) * (hz / hx) * (c(:, 1:end-1) - c(:, 2:end));

    Fx = Fx_adv + Fx_diff;
    Fz = Fz_adv + Fz_diff;
    Fz_top = (1 / Ra) * (hx / hz) * (c_top - 2.*c(1, :));
    Fz_bottom = zeros(1, Nx);
    Fx_left = zeros(Nz, 1);
    Fx_right = zeros(Nz, 1);

    % Fz_bottom = (1 / Ra) * (hx / hz) * (zeros(1, Nx) - 2.*c(end, :));
    % Fx_left = (1 / Ra) * (hz / hx) * (zeros(Nz, 1) - 2.*c(:, 1));
    % Fx_right = (1 / Ra) * (hz / hx) * (zeros(Nz, 1) - 2.*c(end, 1));

    Fx = [Fx_left, Fx, Fx_right];
    Fz = [Fz_top; Fz; Fz_bottom];

    dc = (dt / (hx * hz)) * (Fz(1:end-1,:) - Fz(2:end, :) + Fx(:, 1:end-1) - Fx(:, 2:end));
    c = c + dc;

    % plot concentration
    if (mod(i, drawdelay) == 0)
    subplot(1,2,1);
    surf(X,Z,-ones(size(Z)), c, 'edgecolor', 'none', 'facecolor', 'interp');
    set(gca, 'YDir', 'reverse');
    view([0 0 1]);
    % colormap(hot);
    % caxis([0, 1]);
    colorbar;
    title(sprintf('Time = %g, ct = %g', t, sum(c(:))));
    hold on;
    % draw appoximate velocities on perm field
    scale = 0.1;
    %figure; h = surf(X,Y,zeros(size(X)), log(perm),'edgecolor','none'); view([0, 0, 1]);
    %axis equal square;
    %hold on;
    % quiver(X,Z,scale*(ux(:, 1:end-1)+ux(:, 2:end))/2, scale*(uz(1:end-1, :) + uz(2:end, :))/2, 'Autoscale','off', 'color', [0,0,0]);
    quiver(X,Z,scale*(ux(:, 1:end-1)+ux(:, 2:end))/2, scale*(uz(1:end-1, :) + uz(2:end, :))/2, 'color', [0,0,0]);
    axis equal;
    hold off;
    subplot(1,2,2);
    surf(X,Z,-ones(size(Z)), p, 'edgecolor', 'none', 'facecolor', 'interp');
    axis equal;
    set(gca, 'YDir', 'reverse');
    view([0 0 1]);
    % colormap(hot);
    % caxis([0, 1]);
    colorbar;
    drawnow;
    hold off;
    end

    i = i + 1;
end

end
