function [A, B] = cfd1g(N, dx, alpha, a, b)
    du = alpha * ones(N, 1);
    dd = ones(N, 1);
    dl = du;
    A = spdiags([dl, dd, du], [1 0 -1], N, N)';
    A(1, end) = alpha;
    A(end, 1) = alpha;

    bll = -(b / (4*dx)) * ones(N, 1);
    bl = -(a / (2*dx)) * ones(N, 1);
    bu = (a / (2*dx)) * ones(N, 1);
    buu = (b / (4*dx)) * ones(N, 1);
    B = spdiags([bll, bl, bu, buu], [2 1 -1 -2], N, N)';
    B(1, end) = bl(1);
    B(1, end-1) = bll(1);
    B(end, 1) = bu(1);
    B(end, 2) = buu(1);
end
