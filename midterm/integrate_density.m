function [ cm, ta, c_total, t_draw ]  = integrate_density( alpha, Nx, Nz, Ra, dt, tf, drawdelay, t_wanted )
%INTEGRATE_DENSITY Solves 2D linear scalar transport equation for density driven flow
%using finite volumes.
%  alpha - aspect ratio (length / height)
%  Nx - number of grid cells in x
%  Ny - number of grid cells in x
%  Ra - Rayleigh number
%  dt - time step size
%  tf - max time
%  drawdelay (optinonal) - how many steps to do before drawing concentration.
%  t_wanted (optional) - time to save concentrations at.
% 
%  Returns 
%       Nz x Nx x (# times) matrix cm containing cell concentrations.
%       1 x (# times) matrix ta containing actual times data is saved at.
if (nargin == 6)
    drawdelay = 1000;
    t_wanted = [];
end

if (nargin == 7)
    t_wanted = [];
end

if (isempty(t_wanted))
    cm = [];
    ta = [];
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
j = 1;
c = zeros(Nz, Nx);

% persistent figure
h = figure;
set(h, 'units', 'inches', 'position', [1 1 16 6])
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
while (t < tf)
    t = t + dt;

    c_z_avg = hx*0.5*(c(1:end-1, :) + c(2:end, :));
    c_z = [zeros(1, Nx); c_z_avg; zeros(1, Nx)];
    b = (c_z(1:end-1,:) - c_z(2:end,:)); 

    bv = reshape(b, [], 1);

    % solve for pressure
    clear pv;
    pv(perm,:) = L'\(D\(L\(bv(perm,:))));

    p = reshape(pv, Nz, Nx);

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
    Fz_top = (1 / Ra) * (hx / hz) * (2.0*c_top - 2.*c(1, :));
    Fz_bottom = zeros(1, Nx);
    Fx_left = zeros(Nz, 1);
    Fx_right = zeros(Nz, 1);

    Fx = [Fx_left, Fx, Fx_right];
    Fz = [Fz_top; Fz; Fz_bottom];

    dc = (dt / (hx * hz)) * (Fz(1:end-1,:) - Fz(2:end, :) + Fx(:, 1:end-1) - Fx(:, 2:end));
    c = c + dc;

    % plot concentration, velocity, pressure
    if (mod(i, drawdelay) == 0)
        subplot(1,3,1);
        surf(X,Z,-ones(size(Z)), c, 'edgecolor', 'none', 'facecolor', 'interp');
        title(sprintf('Concentration at Time = %g', t)); 
        set(gca, 'YDir', 'reverse');
        view([0 0 1]);
        % colormap(hot);
        colorbar;
        hold on;
        % draw appoximate velocities on concentration field
        scale = 0.1;
        % quiver(X,Z,scale*(ux(:, 1:end-1)+ux(:, 2:end))/2, scale*(uz(1:end-1, :) + uz(2:end, :))/2, 'Autoscale','off', 'color', [0,0,0]);
        quiver(X,Z,scale*(ux(:, 1:end-1)+ux(:, 2:end))/2, scale*(uz(1:end-1, :) + uz(2:end, :))/2, 'color', [0,0,0]);
        axis equal;
        xlim([0, 1]);
        ylim([0, alpha]);
        caxis([0, 1]);
        hold off;

        subplot(1,3,2);
        surf(X,Z,-ones(size(Z)), p, 'edgecolor', 'none', 'facecolor', 'interp');
        title('Pressure');
        axis equal;
        set(gca, 'YDir', 'reverse');
        view([0 0 1]);
        colorbar;
        xlim([0, 1]);
        ylim([0, alpha]);
        
        subplot(1,3,3);
        quiver(X,Z,scale*(ux(:, 1:end-1)+ux(:, 2:end))/2, scale*(uz(1:end-1, :) + uz(2:end, :))/2, 'color', [0,0,0]);
        title('Velocity');
        axis equal;
        set(gca, 'YDir', 'reverse');
        view([0 0 1]);
        xlim([0, 1]);
        ylim([0, alpha]);

        c_total(i/drawdelay) = hx*hz*sum(c(:));
        t_draw(i/drawdelay) = t;

        drawnow;
        fname = sprintf('figs/Ra_%d_%d.png', Ra, i);
        print(fname, '-dpng');
    end

    if (~isempty(t_wanted))
        if (t > t_wanted(1))
            ta(j) = t;
            cm(:, :, j) = c(:, :);
            t_wanted = t_wanted(2:end);
            j = j + 1;
        end
    end
    i = i + 1;

end

end
