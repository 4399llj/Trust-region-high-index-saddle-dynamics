function r = paper_solve(m, x, o)
    % Algorithm 2.1: fixed-base pair acceptance and independent radii.
    defaults = struct('method', 'CP', 'maxiter', 1e5, 'tol', 1e-5, 'delta', 1e-3, ...
                     'cap', Inf, 'eigtol', 1e-3, 'exactSpectrum', false, ...
                     'tcgmax', 50, 'tcgtol', 1e-3, 'dt', 1e-3, ...
                     'bbmin', 1e-10, 'bbmax', 1e-1);
    f = fieldnames(defaults);
    for i = 1:numel(f)
        if ~isfield(o, f{i})
            o.(f{i}) = defaults.(f{i});
        end
    end
    delta = [o.delta, o.delta];
    V = [];
    if isfield(o, 'V')
        V = o.V;
    end
    n = numel(x);
    t = tic;
    residual = zeros(o.maxiter + 1, 1);
    errors = residual;
    eigres = residual;
    rejects = 0;
    oldx = [];
    oldf = [];
    beta = o.dt;
    step_history = zeros(o.maxiter, 1);
    for iter = 0:o.maxiter
        if mod(iter, 1000) == 0 && iter > 0
            fprintf('%s iteration %d residual %.4g\n', o.method, iter, norm(m.grad(x)));
        end
        g = m.grad(x);
        residual(iter + 1) = norm(g);
        errors(iter + 1) = norm(x - m.exact);
        if norm(g) < o.tol || iter == o.maxiter
            break
        end
        A = m.hessian(x);
        if o.exactSpectrum || strcmp(o.method, 'N')
            [U, L] = eig(full(A), 'vector');
            [lam, j] = sort(L);
            U = U(:, j);
            V = U(:, 1:m.k);
            er = norm(A * V - V * diag(lam(1:m.k)), 'fro');
        else
            if isfield(m, 'precondition')
                pre = m.precondition;
            else
                pre = @(z) z;
            end
            [V, ~, er] = paper_eigenspace(A, m.k, V, o.eigtol, pre);
        end
        eigres(iter + 1) = er;
        proj = @(z) z - V * (V' * z);
        g1 = V * (V' * g);
        g2 = proj(g);
        if any(strcmp(o.method, {'HiSD', 'BB', 'AHiSD'}))
            force = g - 2 * g1;
            if strcmp(o.method, 'BB') && ~isempty(oldx)
                sx = x - oldx;
                y = force - oldf;
                beta = min(o.bbmax, max(o.bbmin, abs(sx' * y) / max(y' * y, realmin)));
            end
            step = -beta * force;
            step_history(iter + 1) = beta;
            if strcmp(o.method, 'AHiSD') && ~isempty(oldx)
                momentum = 0.5;
                if isfield(o, 'momentum')
                    momentum = o.momentum;
                end
                step = step + momentum * (x - oldx);
            end
            oldx = x;
            oldf = force;
            x = x + step;
            continue
        end
        accepted = false;
        for trial = 1:100
            if strcmp(o.method, 'N')
                p1 = V * trustdiag(-lam(1:m.k), V' * g, delta(1));
                W = U(:, m.k + 1:end);
                p2 = W * trustdiag(lam(m.k + 1:end), -W' * g, delta(2));
            else
                p1 = cauchy(g1, -A, delta(1));
                if strcmp(o.method, 'TCG')
                    p2 = tcg(g2, A, proj, delta(2), o.tcgmax, o.tcgtol);
                else
                    p2 = cauchy(-g2, A, delta(2));
                end
            end
            pred = [g' * p1 + p1' * A * p1 / 2, -g' * p2 - p2' * A * p2 / 2];
            actual = [m.diff(x, p1), -m.diff(x, p2)];
            rho = zeros(1, 2);
            for b = 1:2
                if pred(b) == 0
                    rho(b) = 0.5;
                elseif pred(b) > 0
                    rho(b) = actual(b) / pred(b);
                else
                    rho(b) = -Inf;
                end
            end
            if all(rho >= 0) && all(isfinite(rho))
                x = x + p1 + p2;
                delta(rho > 0.75) = 2 * delta(rho > 0.75);
                delta(rho < 0.25) = delta(rho < 0.25) / 4;
                delta = min(delta, o.cap);
                accepted = true;
                break
            end
            delta(rho < 0 | ~isfinite(rho)) = delta(rho < 0 | ~isfinite(rho)) / 4;
            rejects = rejects + 1;
        end
        if ~accepted
            error('paper:rejection', 'Trial pair failed after 100 retries');
        end
    end
    r = struct('x', x, 'iterations', iter, ...
               'residual', residual(1:iter + 1), 'error', errors(1:iter + 1), ...
               'eigres', eigres(1:iter), 'seconds', toc(t), ...
               'converged', norm(m.grad(x)) < o.tol, ...
               'rejections', rejects, 'options', o);
    r.bb_steps = step_history(1:iter);
end

function p = cauchy(d, A, r)
    q = d' * d;
    if q == 0
        p = d;
        return
    end
    c = d' * A * d;
    a = r / sqrt(q);
    if c > 0
        a = min(a, q / c);
    end
    p = a * d;
end

function p = trustdiag(a, b, r)
    % Exact diagonal quadratic trust-region minimizer, including hard case.
    lo = max(0, -min(a));
    tiny = 1e-12 * max(1, max(abs(a)));
    den = a + lo;
    idx = abs(den) < tiny;
    p = zeros(size(b));
    p(~idx) = b(~idx) ./ den(~idx);
    if all(abs(b(idx)) < tiny) && norm(p) <= r
        if lo > 0
            j = find(idx, 1);
            p(j) = sqrt(max(0, r * r - p' * p));
        end
        return
    end
    hi = max(1, lo + 1);
    while norm(b ./ (a + hi)) > r
        hi = 2 * hi + 1;
    end
    for j = 1:100
        mid = (lo + hi) / 2;
        if norm(b ./ (a + mid)) > r
            lo = mid;
        else
            hi = mid;
        end
    end
    p = b ./ (a + hi);
end

function p = tcg(g, A, proj, r, maxit, tol)
    p = zeros(size(g));
    z = g;
    d = -z;
    rr = z' * z;
    rr0 = rr;
    if rr == 0
        return
    end
    for j = 1:maxit
        Ad = proj(A * d);
        c = d' * Ad;
        tau = (-p' * d + sqrt((p' * d)^2 + (d' * d) * (r * r - p' * p))) / (d' * d);
        if c <= 0
            p = p + tau * d;
            return
        end
        a = rr / c;
        if norm(p + a * d) >= r
            p = p + tau * d;
            return
        end
        p = p + a * d;
        z = proj(z + a * Ad);
        rn = z' * z;
        if sqrt(rn) <= tol * sqrt(rr0)
            return
        end
        d = -z + (rn / rr) * d;
        rr = rn;
    end
end
