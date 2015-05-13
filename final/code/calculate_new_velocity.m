function [un, psi_tilde_hat] = calculate_new_velocity(cn, uo, ubg, Lx, Ly, R)
    % velocity stored as u(xpos, ypos, component).
    % u(:,:,1) are ux for entire grid, u(:,:,2) are uy for entire grid
    Nx = size(uo, 2);
    Ny = size(uo, 1);
    kx_odd= (2*pi / Lx) * [0:Nx/2-1, 0,-Nx/2+1:-1];
    ky_odd = (2*pi / Ly) * [0:Ny/2-1, 0, -Ny/2+1:-1];
    kx_even = (2*pi / Lx) * [0:Nx/2, -Nx/2+1:-1];
    ky_even = (2*pi / Ly) * [0:Ny/2, -Ny/2+1:-1];
    [mkx_odd, mky_odd] = meshgrid(kx_odd, ky_odd);
    [mkx_even, mky_even] = meshgrid(kx_even, ky_even);
    c_hat = fft2(cn);
    dcdx = real(ifft2(1i*mkx_odd.*c_hat));
    dcdy = real(ifft2(1i*mky_odd.*c_hat));
    % figure(3);
    % subplot(2,1,1);
    % surf(dcdx, 'edgecolor', 'none');
    % colorbar;
    % view([0 0 1]);
    % subplot(2,1,2);
    % surf(dcdy, 'edgecolor', 'none');
    % colorbar;
    % view([0 0 1]);
    % drawnow;
    % max(dcdx(:))
    % max(dcdy(:))
    omega = R*(dcdx.*(uo(:,:,2) + ubg(:,:,2)) - dcdy.*(uo(:,:,1) + ubg(:,:,1)));
    omega_hat = fft2(omega);
    imksq = -1 ./ (mkx_even.^2 + mky_even.^2);
    imksq(1,1) = 1; % make sure 'free' average value doesn't cause issues
    psi_tilde_hat = imksq .* -omega_hat;
    un = zeros(size(uo));
    un(:,:,1) = real(ifft2(1i*mky_odd.*psi_tilde_hat));
    un(:,:,2) = -real(ifft2(1i*mkx_odd.*psi_tilde_hat));
    % max(un(:))
end
