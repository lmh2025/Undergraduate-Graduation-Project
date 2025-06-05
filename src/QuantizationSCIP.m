function X = QuantizationSCIP(A, Xideal)
    [m, n] = size(A);
    H = spdiags(1 ./ Xideal(:).^2, 0, n, n); f = -2 ./ Xideal(:);
    lhs = zeros(m, 1); rhs = zeros(m, 1);
    lb = ones(n, 1); ub = inf(n, 1);
    vtype = repmat('I', 1, n);

    X = scip(H, f, A, lhs, rhs, lb, ub, vtype);
    X = round(X);
end