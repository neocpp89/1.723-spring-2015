clear all;
close all;

L = 1; N = 256;
x = 0:L/N:L; x = x(1:N); y = x;
dx = x(2) - x(1);
[xx,yy] = meshgrid(x,y);
rnd = (-1/N)+(2/N)*rand(1,N);
yr = repmat((y+rnd)',1,N); % Perturbed y coordinate
w = [20 21 30 31 40]; % Fourier modes in the perturbed front
p = 0*yy;
for iw = 1:length(w)
p = p + 0.004*(rand-0.5)*cos(w(iw)*pi*yr+rand*iw);
end
c = 0.5*(erfc(100*(p+xx-0.75*L)) - erfc(100*(p+xx-0.25*L)));

figure(1);
colormap('hot');
surf(xx,yy,-ones(size(xx)), c, 'edgecolor', 'none', 'facecolor', 'interp');
caxis([0 1]);
view([0 0 1]);
axis equal square;

Pe = 2500;
R = 1;

t = 0;
tf = 4;
u = zeros([size(c) 2]);
ubg = zeros(size(u));
ubg(:, :, 1) = 1;

step = 0;

figure(2);
while t < tf
    umag = sqrt(squeeze(u(:,:,1)).^2 + squeeze(u(:,:,2)).^2);
    dt_adv = dx / (max(1e-8, max(umag(:))));
    dt_diff = Pe * dx * dx;
    dt = 0.1 * min(dt_adv, dt_diff);

    c = calculate_new_concentration(c, u, ubg, L, L, Pe, dt);
    [u, psi_tilde_hat] = calculate_new_velocity(c, u, ubg, L, L, R);

    if (mod(step, 50) == 0)
        subplot(2,1,1);
        surf(xx,yy,-ones(size(xx)), c, 'edgecolor', 'none', 'facecolor', 'interp');
        title(sprintf('t=%g', t));
        % caxis([0 1]);
        colorbar;
        view([0 0 1]);
        axis equal square;
        drawnow;

        subplot(2,1,2);
        ux(:,:) = u(:,:,1) + ubg(:,:,1); 
        uy(:,:) = u(:,:,2) + ubg(:,:,2); 
        ux(:,:) = u(:,:,1); 
        uy(:,:) = u(:,:,2); 
        quiver(xx,yy,ux,uy);
        axis equal square;
        drawnow;
    end

    t = t + dt;
    step = step + 1;
end

% var(c(:));
