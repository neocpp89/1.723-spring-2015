function [ XYm, Cm, ta, p ] = integrate_scalar_tracer_equation_2d( Nx, Ny, var_lnk, Pe, theta, dt, tout )
%INTEGRATE_SCALAR_TRACER_EQUATION_2D Solves 2D linear scalar transport equation
%using finite volumes.
%  Nx - number of grid cells in x
%  Ny - number of grid cells in x
%  Pe - Peclet number
%  theta - trapezoidal integration parameter (0 = forward euler, 1 = backward euler) 
%  dt - time step size
%  tout - list of times to save solution vector
% 
%  Returns Nx x Ny x (# times) matrix xm containing cell point locations,
%          Nx x Ny x (# times) matrix cm containing cell solutions at different times,
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

lambda_bar_interior_x = (hy / hx) * 2 ./ (kinv(:,2:Ny) + kinv(:,1:Ny-1));
lambda_bar_x = zeros(Nx,Ny+1);
lambda_bar_x(:,2:end-1) = lambda_bar_interior_x(:,:);
lambda_bar_interior_y = (hx / hy) * 2 ./ (kinv(2:Nx,:) + kinv(1:Nx-1,:));
lambda_bar_y = zeros(Nx+1,Ny);
lambda_bar_y(2:end-1,:) = lambda_bar_interior_y(:,:);


% Set BCs (default is no-flow)
b = zeros(Nx,Ny);
b(1) = 1;

lambda_bar_x(end,end) = (hy / hx) * 2 * perm(end,end);
lambda_bar_y(end,end) = (hx / hy) * 2 * perm(end,end);

% make diagonals
T_x_left = reshape(lambda_bar_x(:, 1:end-1), [], 1);
T_x_right = reshape(lambda_bar_x(:, 2:end), [], 1);
T_y_bottom = reshape(lambda_bar_y(1:end-1, :), [], 1);
T_y_top = reshape(lambda_bar_y(2:end, :), [], 1);

sumT = T_x_right + T_x_left + T_y_bottom + T_y_top;
diags = [-T_x_left, -T_y_bottom, sumT, -T_y_top, -T_x_right];

% assemble matrix
A = spdiags(diags, [Nx;1;0;-1;-Nx], Nx*Ny, Nx*Ny)';
bv = reshape(b, [], 1);

% solve for pressure
pv = A \ bv;
p = reshape(pv, Nx, Ny);

figure; surf(X,Y,log(perm),'edgecolor','none'); view([0, 0, 1]);
figure; surf(X,Y,p,'edgecolor','none'); view([0, 0, 1]);
figure; plot(diag(p) - p(end/2,end/2));

ux = lambda_bar_x(:, 2:end-1).*(p(:, 1:end-1) - p(:,2:end));
uy = lambda_bar_y(2:end-1, :).*(p(1:end-1,:) - p(2:end,:));
zx = zeros(size(ux, 1), 1);
zy = zeros(1, size(uy, 2));
ux = [zx, ux, zx];
uy = [zy; uy; zy];
ux(1,1) = -1/2;
uy(1,1) = -1/2;
ux(end,end) = 1/2;
uy(end,end) = 1/2;

scale = 0.1;
figure; h = surf(X,Y,zeros(size(X)), log(perm),'edgecolor','none'); view([0, 0, 1]);
axis equal square;
hold on;
quiver(X,Y,scale*(ux(:, 1:end-1)+ux(:, 2:end))/2, scale*(uy(1:end-1, :) + uy(2:end, :))/2, 'Autoscale','off', 'color', [0,0,0]);
hold off;

e = ones(N, 1);
z = zeros(N, 1);
alpha = (1.0/ (Pe * h));
diags = (dt / h) * ([-e e z] - alpha * [e -2*e e]);

tout = sort(tout);

% initial condition
c = zeros(N, 1);

c_tilde = zeros(N, 1);

% left BC
diags(1, :) = (dt / h) * [0 (1+alpha) -alpha];
c_tilde(1) = dt/h * (1);

% right BC
diags(end, :) = (dt / h) * [-(1+alpha) (1+alpha) 0];

A = spdiags(diags, [1;0;-1], N, N)';
I = speye(N, N);

LA = I+A*theta;
RA = I-A*(1.0-theta);

xm = kron(ones(numel(tout), 1), x)';
cm = zeros(N, numel(tout));
t = 0;
i = 1;
ta = zeros(numel(tout), 1);
while i <= numel(tout)
    if (theta == 0)
        c = RA*c + c_tilde;
    else
        c = LA \ (RA*c + c_tilde);
    end
    % plot(x,c);
    % ylim([0, max(1,max(c))]);
    % drawnow;
    t = t + dt;
    if (t >= tout(i))
        ta(i) = t;
        cm(:, i) = c;
        i = i + 1;
    end
end


end
