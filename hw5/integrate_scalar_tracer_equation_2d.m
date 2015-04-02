function [ Cm, ta, p ] = integrate_scalar_tracer_equation_2d( Nx, Ny, var_lnk, Pe, theta, dt, tf )
%INTEGRATE_SCALAR_TRACER_EQUATION_2D Solves 2D linear scalar transport equation
%using finite volumes.
%  Nx - number of grid cells in x
%  Ny - number of grid cells in x
%  Pe - Peclet number
%  theta - trapezoidal integration parameter (0 = forward euler, 1 = backward euler) 
%  dt - time step size
%  tf - max time
% 
%  Returns Nx x Ny x (# times) matrix xm containing cell point locations,
%          Nx x Ny x (# times) matrix cm containing cell concentrations at different times,
%          1 x (# times) matrix ta containing actual times data is saved at.
%          Nx x Ny matrix containing pressure field.

Lx = 1;
Ly = 1;
hx = Lx / Nx;
hy = Ly / Ny;
x = linspace(0+hx/2.0, Lx-hx/2.0, Nx);
y = linspace(0+hy/2.0, Ly-hy/2.0, Ny);
[X,Y] = meshgrid(x,y);
X = X';
Y = Y';

corr_lenx = 4*hx;
corr_leny = 4*hy;

[perm,var_lnk_actual]= random_perm(var_lnk,corr_lenx,corr_leny,Nx,Ny,Lx,Ly);
perm = perm';
kinv = 1./perm;

lambda_bar_x = zeros(Nx+1,Ny);
lambda_bar_x(2:end-1,:) = (hy / hx) * 2 ./ (kinv(2:Nx,:) + kinv(1:Nx-1,:));
lambda_bar_y = zeros(Nx,Ny+1);
lambda_bar_y(:,2:end-1) = (hx / hy) * 2 ./ (kinv(:,2:Ny) + kinv(:,1:Ny-1));


% Set BCs (default is no-flow)
b = zeros(Nx,Ny);
b(1,1) = 1;

lambda_bar_x(end,end) = (hy / hx) * 2 * perm(end,end);
lambda_bar_y(end,end) = (hx / hy) * 2 * perm(end,end);

% make diagonals
T_x_left = reshape(lambda_bar_x(1:end-1, :)', [], 1);
T_x_right = reshape(lambda_bar_x(2:end, :)', [], 1);
T_y_bottom = reshape(lambda_bar_y(:, 1:end-1)', [], 1);
T_y_top = reshape(lambda_bar_y(:, 2:end)', [], 1);

sumT = T_x_right + T_x_left + T_y_bottom + T_y_top;
diags = [-T_x_left, -T_y_bottom, sumT, -T_y_top, -T_x_right];

% assemble matrix
A = spdiags(diags, [Ny;1;0;-1;-Ny], Nx*Ny, Nx*Ny)';
bv = reshape(b, [], 1);

% solve for pressure
pv = A \ bv;
p = reshape(pv, Nx, Ny);

% show plots and save pressure
figure; surf(X,Y,log(perm),'edgecolor','none'); view([0, 0, 1]);
figure; surf(X,Y,p,'edgecolor','none'); view([0, 0, 1]);
h = figure; plot(x, diag(p) - p(end/2,end/2));
set(h, 'units', 'inches', 'position', [1 1 3 3])
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
fname = sprintf('figs/pressure_%d_%d.png', 10*var_lnk, Nx);
title(sprintf('Pressure along diagonal\n(N = %d)', Nx));
xlabel('Coordinate');
ylabel('Pressure');
print(fname, '-dpng');

% calculate integrated velocities
ux_int = lambda_bar_x(2:end-1, :).*(p(1:end-1, :) - p(2:end, :));
uy_int = lambda_bar_y(:, 2:end-1).*(p(:, 1:end-1) - p(:, 2:end));
zx = zeros(1, size(ux_int, 2));
zy = zeros(size(uy_int, 1), 1);
ux = [zx; ux_int; zx];
uy = [zy, uy_int, zy];
ux(1,1) = 1/2;
uy(1,1) = 1/2;
ux(end,end) = 1/2;
uy(end,end) = 1/2;

% draw appoximate velocities on perm field
scale = 0.1;
figure; h = surf(X,Y,zeros(size(X)), log(perm),'edgecolor','none'); view([0, 0, 1]);
axis equal square;
hold on;
quiver(X,Y,scale*(ux(1:end-1,:)+ux(2:end, :))/2, scale*(uy(:, 1:end-1) + uy(:, 2:end))/2, 'Autoscale','off', 'color', [0,0,0]);
hold off;

t = 0;
i = 1;
c = zeros(Nx, Ny);
while (t < tf)
    t = t + dt;
    ta(i) = t;
    Cm(:,:, i) = c;
    c(1,1) = 1;

    Fx_adv = ux_int.*(ux_int > 0).*c(1:end-1, :) + ux_int.*(ux_int < 0).*c(2:end, :);
    zFx = zeros(1, size(Fx_adv, 2));
    Fx_adv = [zFx; Fx_adv; zFx];

    Fy_adv = uy_int.*(uy_int > 0).*c(:, 1:end-1) + uy_int.*(uy_int < 0).*c(:, 2:end);
    zFy = zeros(size(Fy_adv, 1), 1);
    Fy_adv = [zFy, Fy_adv, zFy];

    Fx_diff = (1 / Pe) * (hy / hx) * (c(1:end-1, :) - c(2:end, :));
    Fx_diff = [zFx; Fx_diff; zFx];

    Fy_diff = (1 / Pe) * (hx / hy) * (c(:, 1:end-1) - c(:, 2:end));
    Fy_diff = [zFy, Fy_diff, zFy];

    Fx = Fx_adv + Fx_diff;
    Fy = Fy_adv + Fy_diff;

    Fx(1,1) = 1/2;
    Fy(1,1) = 1/2;
    Fx(end,end) = c(end,end)*1/2;
    Fy(end,end) = c(end,end)*1/2;

    dc = (dt / (hx * hy)) * (Fx(1:end-1,:) - Fx(2:end, :) + Fy(:, 1:end-1) - Fy(:, 2:end));
    c = c + dc;

    % plot concentration
    % surf(X,Y,c, 'edgecolor', 'none');
    % view([0 0 1]);
    % colormap(hot);
    % drawnow;

    i = i + 1;
end

end
