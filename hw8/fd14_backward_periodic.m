function [D] = fd14_backward_periodic(N, dx)
    z = zeros(1, N);
    r = -[-25/12 4 -3  4/3 -1/4] ./ dx;
    c = z;
    z(1:numel(r)) = r;
    c(1) = r(1);
    cc = numel(c) - numel(r(2:end)) + 1;
    c(cc:end) = r(end:-1:2);
    D = sparse(toeplitz(z,c));
end
