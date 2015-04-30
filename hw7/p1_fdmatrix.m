function [D] = p1_fdmatrix(N)
    ll = local_fdcoeffs([0 -1 -2 -3], 1)
    l = local_fdcoeffs([1 0 -1 -2], 1)
    c = local_fdcoeffs([2 1 -1 -2], 1)
    r = local_fdcoeffs([2 1 0 -1], 1)
    rr = local_fdcoeffs([3 2 1 0], 1)
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
    A
    fdc = A \ b;
end
