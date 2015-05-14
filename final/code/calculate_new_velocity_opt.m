function [u_tilde_new, psi_tilde_hat] = calculate_new_velocity_opt(c_new, u_tilde_old, u_background, R, mkx_odd, mky_odd, mksq)
    % Velocity is stored as u(xpos, ypos, component).
    % i.e. u(:,:,1) are ux for entire grid, u(:,:,2) are uy for entire grid
    c_hat = fft2(c_new);
    dcdx = real(ifft2(1i*mkx_odd.*c_hat));
    dcdy = real(ifft2(1i*mky_odd.*c_hat));

    omega = R*(dcdx.*(u_tilde_old(:,:,2) + u_background(:,:,2)) - dcdy.*(u_tilde_old(:,:,1) + u_background(:,:,1)));
    omega_hat = fft2(omega);

    psi_tilde_hat = -omega_hat ./ mksq;
    u_tilde_new = zeros(size(u_tilde_old));
    u_tilde_new(:,:,1) = real(ifft2(1i*mky_odd.*psi_tilde_hat));
    u_tilde_new(:,:,2) = -real(ifft2(1i*mkx_odd.*psi_tilde_hat));
end
