function [V, lambda, res] = paper_eigenspace(A, k, V, tol, precondition)
    % Warm-start block locally optimal Rayleigh-Ritz (LOBPCG, no toolbox).
    n = size(A, 1);
    if nargin < 5
        precondition = @(r) r;
    end
    if k == 0
        V = zeros(n, 0);
        lambda = [];
        res = 0;
        return
    end
    if n <= 2 * k + 1
        [U, L] = eig(full(A), 'vector');
        [lambda, j] = sort(L);
        V = U(:, j(1:k));
        lambda = lambda(1:k);
        res = norm(A * V - V * diag(lambda), 'fro');
        return
    end
    if isempty(V)
        V = randn(n, k);
    end
    [V, ~] = qr(V, 0);
    P = [];
    for it = 1:1000
        [U, L] = eig(full(V' * A * V), 'vector');
        [lambda, j] = sort(L);
        V = V * U(:, j);
        R = A * V - V * diag(lambda);
        res = norm(R, 'fro');
        if res <= tol
            return
        end
        Z = [V, precondition(R), P];
        [U, S, ~] = svd(Z, 'econ');
        sv = diag(S);
        Z = U(:, sv > max(sv) * 1e-12);
        [C, L] = eig(full(Z' * A * Z), 'vector');
        [~, j] = sort(L);
        T = Z * C(:, j(1:k));
        P = T - V * (V' * T);
        V = T;
    end
    error('paper:eigensolver', 'LOBPCG residual %.3g exceeds tolerance %.3g', res, tol);
end
