function [D] = p1_fdmatrix(N)
% returns the sparse NxN matrix for differentiating problem 1.
    ll = local_fdcoeffs([0 1 2 3], 1);
    l = local_fdcoeffs([-1 0 1 2], 1);
    c = local_fdcoeffs([-2 -1 1 2], 1);
    r = local_fdcoeffs([-2 -1 0 1], 1);
    rr = local_fdcoeffs([-3 -2 -1 0], 1);
    cols = zeros(N, 5);
    C = kron(ones(N-4,1),[c(1) c(2) 0 c(3) c(4)]);
    cols(3:end-2, :) = C;
    D = spdiags(cols, [2 1 0 -1 -2], N, N)';
    D(1, 1:4) = ll;
    D(2, 1:4) = l;
    D(end-1, end-3:end) = r;
    D(end, end-3:end) = rr;
    % dont forget to divide by h (= 1/(N-1))
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
