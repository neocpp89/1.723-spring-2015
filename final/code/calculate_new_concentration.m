function [cn] = calculate_new_concentration(co, uo, ubg, Lx, Ly, Pe, dt)
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
    mksq = -(mkx_even.^2 + mky_even.^2);

    Y1 = f(co, uo, ubg, Pe, mkx_odd, mky_odd, mksq);
    Y2 = f(co + dt*0.5*Y1, uo, ubg, Pe, mkx_odd, mky_odd, mksq);
    Y3 = f(co + dt*0.5*Y2, uo, ubg, Pe, mkx_odd, mky_odd, mksq);
    Y4 = f(co + dt*Y3, uo, ubg, Pe, mkx_odd, mky_odd, mksq);
    cn = co + (1/6)*dt*(Y1 + 2*Y2 + 2*Y3 + Y4);
end

function [Y] = f(co, uo, ubg, Pe, mkx_odd, mky_odd, mksq)
    uxc = squeeze((uo(:,:,1))+ubg(:,:,1)).*co;
    uyc = squeeze((uo(:,:,2))+ubg(:,:,2)).*co;
    uxc_hat = fft2(uxc);
    uyc_hat = fft2(uyc);
    duxcdx = real(ifft2(1i*mkx_odd.*uxc_hat));
    duycdy = real(ifft2(1i*mky_odd.*uyc_hat));

    c_hat = fft2(co);
    dcsq = real(ifft2(mksq .* c_hat));
    Y = -(duxcdx + duycdy) + (1 / Pe) * dcsq;
end
