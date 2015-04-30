function [D] = p2_fdmatrix(N)
% returns the sparse NxN matrix for differentiating problem 2.
    cc = local_fdcoeffs([-2 -1 0 1 2], 2);
    r = zeros(1, N);
    r(1) = cc(3);
    r(2) = cc(4);
    r(3) = cc(5);
    r(end-1) = cc(1);
    r(end) = cc(2);
    D = sparse(toeplitz(r));
    % dont forget to divide by h^2 (= 4*pi^2/((N-1)(N-1)))
end

function [c] = taylor_coeffs(di, nterms)
% di is the difference in index from the expansion point
% nterms is how many terms to expand out (length of c)
    c = zeros(nterms, 1);
    c(1) = 1;
    for i=1:(nterms-1)
        c(i+1) = (di ^ i)  / factorial(i);
    end
end

function [fdc] = local_fdcoeffs(vdi, whichderiv)
    nterms = numel(vdi);
    b = zeros(nterms, 1);
    b(whichderiv+1) = 1;
    for i=1:numel(vdi)
        A(:, i) = taylor_coeffs(vdi(i), nterms);
    end
    fdc = A \ b;
end
