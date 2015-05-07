function [D] = fd24_central_periodic(N, dx)
    z = zeros(1, N);

    r = [-1/12 4/3 -5/2 4/3 -1/12] ./ (dx*dx);
    c = z;
    z(1:3) = r(3:5);
    z(end-1:end) = r(1:2);
    D = sparse(toeplitz(z,z));
end
