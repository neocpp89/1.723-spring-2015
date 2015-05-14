clear all;
close all;

% scan over parameter
R = [1, 2, 2.8];

% common parameters
tf = 5;
Pe = 2500;
L = 1;
N = 256;

% grid in real space
x = 0:L/N:L;
x = x(1:N);
y = x;
dx = x(2) - x(1);
[xx,yy] = meshgrid(x,y);

% grid in fourier space
Lx = L;
Ly = L;
Nx = N;
Ny = N;
kx_odd= (2*pi / Lx) * [0:Nx/2-1, 0,-Nx/2+1:-1];
ky_odd = (2*pi / Ly) * [0:Ny/2-1, 0, -Ny/2+1:-1];
kx_even = (2*pi / Lx) * [0:Nx/2, -Nx/2+1:-1];
ky_even = (2*pi / Ly) * [0:Ny/2, -Ny/2+1:-1];
[mkx_odd, mky_odd] = meshgrid(kx_odd, ky_odd);
[mkx_even, mky_even] = meshgrid(kx_even, ky_even);
mksq = -(mkx_even.^2 + mky_even.^2);
mksq(1,1) = 1; % make sure 'free' average value doesn't cause issues

% create initial perturbed concentration distribution
rnd = (-1/N)+(2/N)*rand(1,N);
yr = repmat((y+rnd)',1,N); % Perturbed y coordinate
w = [20 21 30 31 40]; % Fourier modes in the perturbed front
p = 0*yy;
for iw = 1:length(w)
p = p + 0.004*(rand-0.5)*cos(w(iw)*pi*yr+rand*iw);
end
c_initial= 0.5*(erfc(100*(p+xx-0.75*L)) - erfc(100*(p+xx-0.25*L)));

% show initial concentration distribution
h = figure;
set(h, 'units', 'inches', 'position', [1 1 5 5])
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
surf(xx,yy,-ones(size(xx)), c_initial, 'edgecolor', 'none', 'facecolor', 'interp');
colorbar;
caxis([0 1]);
view([0 0 1]);
axis equal square;
title('Initial Concentration Field');
xlabel('X');
ylabel('Y');
print('../report/initial_conc.png', '-dpng');

% set backround velocity as (1, 0) everywhere.
u_background = zeros([size(c_initial) 2]);
u_background(:, :, 1) = 1;


tlist = [0.25, 0.5, 1, 2, 4];

% I don't know if this is always available, so it is off by default.
% parfor i=1:numel(R)
for i=1:numel(R)
    j = 1;
    t = 0;
    step = 0;
    c = c_initial;
    u_tilde = zeros(size(u_background));
    figure;
    c_hat = fft2(c);
    dcdx = real(ifft2(1i * mkx_odd .* c_hat));
    dcdy = real(ifft2(1i * mky_odd .* c_hat));
    epsilon_c_saved(1) = mean((1 / Pe) * (dcdx(:).^2 + dcdy(:).^2));
    varc_saved(1) = var(c(:));
    num_fingers_saved(1) = 0;
    t_saved(1) = 0;
    while t < tf
        % set adaptive time step
        umag = sqrt(squeeze(u_tilde(:,:,1) + u_background(:,:,1)).^2 + squeeze(u_tilde(:,:,2) + u_background(:,:,2)).^2);
        dt_adv = dx / (max(1e-8, max(umag(:)))); % avoid divide by zero
        dt_diff = Pe * dx * dx;
        dt = 0.25 * min(dt_adv, dt_diff);

        c = calculate_new_concentration_opt(c, u_tilde, u_background, Pe, dt, mkx_odd, mky_odd, mksq);
        [u_tilde, psi_tilde_hat] = calculate_new_velocity_opt(c, u_tilde, u_background, R(i), mkx_odd, mky_odd, mksq);

        % uncomment to get concentration and periodic velocity plots every few steps.
        %{
        if (mod(step, 50) == 0)
            subplot(1,2,1);
            surf(xx,yy,-ones(size(xx)), c, 'edgecolor', 'none', 'facecolor', 'interp');
            title(sprintf('t=%g', t));
            caxis([0 1]);
            colorbar;
            view([0 0 1]);
            xlim([0, L]);
            ylim([0, L]);
            axis equal square;
            drawnow;

            subplot(1,2,2);
            ux(:,:) = u_tilde(:,:,1); 
            uy(:,:) = u_tilde(:,:,2); 
            quiver(xx,yy,ux,uy);
            xlim([0, L]);
            ylim([0, L]);
            axis equal square;
            drawnow;
        end
        %}

        t = t + dt;
        step = step + 1;

        c_hat = fft2(c);
        dcdx = real(ifft2(1i * mkx_odd .* c_hat));
        dcdy = real(ifft2(1i * mky_odd .* c_hat));
        c_haty = fft(c, [], 1);
        c_haty(1, :) = 0; % remove average value
        powerspec = c_haty.*conj(c_haty);
        powerspec = powerspec(1:end/2, :); % remove aliased values
        [val, idx] = max(powerspec(:));
        [ir, ic] = ind2sub(size(powerspec), idx);
        dominant_k = ky_odd(ir);

        if (step > numel(t_saved))
            epsilon_c_saved = [epsilon_c_saved; zeros(size(epsilon_c_saved))];
            varc_saved = [varc_saved; zeros(size(varc_saved))];
            t_saved = [t_saved; zeros(size(t_saved))];
            num_fingers_saved = [num_fingers_saved; zeros(size(num_fingers_saved))];
        end
        epsilon_c_saved(step) = mean((1 / Pe) * (dcdx(:).^2 + dcdy(:).^2));
        varc_saved(step) = var(c(:));
        t_saved(step) = t;
        num_fingers_saved(step) = dominant_k/(2*pi);

        % show current time every once in a while
        if (mod(step, 50) == 0)
            fprintf('[R=%g] Time: %-5.5f/%-5.5f Current step: %d\n', R(i), t, min(tf,tlist(end)), step);
        end

        % save figures at various times
        if (tlist(j) < t)
            h = figure;
            set(h, 'units', 'inches', 'position', [1 1 4 4])
            set(h, 'PaperUnits','centimeters');
            set(h, 'Units','centimeters');
            pos=get(h,'Position');
            set(h, 'PaperSize', [pos(3) pos(4)]);
            set(h, 'PaperPositionMode', 'manual');
            set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
            surf(xx,yy,-ones(size(xx)), c, 'edgecolor', 'none', 'facecolor', 'interp');
            colorbar;
            caxis([0 1]);
            view([0 0 1]);
            axis equal square;
            title(sprintf('Concentration Field at t=%g\n(R=%g)', t, R(i)));
            xlabel('X');
            ylabel('Y');
            print(sprintf('../report/conc%d_%d.png', 10*R(i), 100*tlist(j)), '-dpng');

            j = j + 1;
            if (j > numel(tlist))
                % That was the last time to save figures, so break.
                % But first, plot statistical data.
                epsilon_c_saved = epsilon_c_saved(1:step);
                t_saved = t_saved(1:step);
                varc_saved = varc_saved(1:step);
                num_fingers_saved = num_fingers_saved(1:step);

                h = figure;
                set(h, 'units', 'inches', 'position', [1 1 9 3])
                set(h, 'PaperUnits','centimeters');
                set(h, 'Units','centimeters');
                pos=get(h,'Position');
                set(h, 'PaperSize', [pos(3) pos(4)]);
                set(h, 'PaperPositionMode', 'manual');
                set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);

                subplot(1,3,1);
                plot(t_saved, varc_saved, 'DisplayName', '\sigma^2(c)');
                legend(gca, 'show', 'location', 'best');
                title(sprintf('Variance\n(R=%g)', R(i)));
                xlabel('Time');
                ylabel('Variance of Concentration');
                subplot(1,3,2);
                plot(t_saved, epsilon_c_saved, 'DisplayName', '\epsilon_c');
                legend(gca, 'show', 'location', 'best');
                title(sprintf('Mixing Rate\n(R=%g)', R(i)));
                xlabel('Time');
                ylabel('Mean Mixing Rate');
                subplot(1,3,3);
                plot(t_saved, num_fingers_saved, 'DisplayName', 'N_f');
                legend(gca, 'show', 'location', 'best');
                title(sprintf('Number of Fingers\n(R=%g)', R(i)));
                xlabel('Time');
                ylabel('Est. Number of Fingers');
                print(sprintf('../report/concstats_%d.png', 10*R(i)), '-dpng');
                break;
            end
        end
    end
end
