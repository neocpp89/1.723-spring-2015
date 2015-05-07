function [A, B] = cfd2g(N, dx, alpha, a, b)
    du = alpha * ones(N, 1);
    dd = ones(N, 1);
    dl = du;
    A = spdiags([dl, dd, du], [1 0 -1], N, N)';
    A(1, end) = alpha;
    A(end, 1) = alpha;

    bll = (b / (4*dx*dx)) * ones(N, 1);
    bl = (a / (dx*dx)) * ones(N, 1);
    bd = -(2*a/(dx*dx) + 2*b/(4*dx*dx)) * ones(N, 1);
    bu = bl;
    buu = bll;
    B = spdiags([bll, bl, bd, bu, buu], [2 1 0 -1 -2], N, N)';
    B(1, end) = bl(1);
    B(1, end-1) = bll(1);
    B(end, 1) = bu(1);
    B(end, 2) = buu(1);
end
