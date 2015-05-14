function [c_new] = calculate_new_concentration_opt(c_old, u_tilde_old, u_background, Pe, dt, mkx_odd, mky_odd, mksq)
    % Implement the RK4 explicit time integration scheme
    Y1 = f(c_old            , u_tilde_old, u_background, Pe, mkx_odd, mky_odd, mksq);
    Y2 = f(c_old + dt*0.5*Y1, u_tilde_old, u_background, Pe, mkx_odd, mky_odd, mksq);
    Y3 = f(c_old + dt*0.5*Y2, u_tilde_old, u_background, Pe, mkx_odd, mky_odd, mksq);
    Y4 = f(c_old + dt*Y3    , u_tilde_old, u_background, Pe, mkx_odd, mky_odd, mksq);
    c_new = c_old + (1/6)*dt*(Y1 + 2*Y2 + 2*Y3 + Y4);
end

function [Y] = f(c_old, u_tilde_old, u_background, Pe, mkx_odd, mky_odd, mksq)
    uxc = squeeze((u_tilde_old(:,:,1))+u_background(:,:,1)).*c_old;
    uyc = squeeze((u_tilde_old(:,:,2))+u_background(:,:,2)).*c_old;
    uxc_hat = fft2(uxc);
    uyc_hat = fft2(uyc);
    duxcdx = real(ifft2(1i*mkx_odd.*uxc_hat));
    duycdy = real(ifft2(1i*mky_odd.*uyc_hat));

    c_hat = fft2(c_old);
    dcsq = real(ifft2(mksq .* c_hat));
    Y = -(duxcdx + duycdy) + (1 / Pe) * dcsq;
end
