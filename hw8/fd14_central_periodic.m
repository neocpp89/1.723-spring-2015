function [D] = fd14_central_periodic(N, dx)
    z = zeros(1, N);
    r = [1/12 -3/4 0 3/4 -1/12] ./ (dx);
    c = z;
    z(1:3) = r(3:5);
    z(end-1:end) = r(1:2);
    c(1:3) = r(3:-1:1);
    c(end-1:end) = r(5:-1:4);
    D = sparse(toeplitz(c,z));
end
