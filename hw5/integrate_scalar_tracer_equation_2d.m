function [ XYm, Cm, ta ] = integrate_scalar_tracer_equation_2d( Nx, Ny, Pe, theta, dt, tout )
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

hx = 1.0 / Nx;
hy = 1.0 / Ny;
x = linspace(0+hx/2.0, 1-hx/2.0, Nx);
y = linspace(0+hy/2.0, 1-hy/2.0, Ny);

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
