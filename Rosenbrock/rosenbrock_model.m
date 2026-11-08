function m = rosenbrock_model(d, k, sneg)
    % Article section 4.1. Analytic sparse Hessian, configurable index.
    s = ones(d, 1);
    s(1:k) = sneg;
    m.energy = @E;
    m.grad = @G;
    m.hessian = @H;
    m.diff = @D;
    m.exact = ones(d, 1);
    m.k = k;

    function e = E(x)
        e = sum(100 * (x(2:end) - x(1:end - 1).^2).^2 ...
                + (1 - x(1:end - 1)).^2) + sum(s .* atan(x - 1).^2);
    end

    function g = G(x)
        r = x(2:end) - x(1:end - 1).^2;
        g = 2 * s .* atan(x - 1) ./ (1 + (x - 1).^2);
        g(1:end - 1) = g(1:end - 1) - 400 * x(1:end - 1) .* r + 2 * (x(1:end - 1) - 1);
        g(2:end) = g(2:end) + 200 * r;
    end

    function A = H(x)
        y = x - 1;
        a = 2 * s .* (1 - 2 * y .* atan(y)) ./ (1 + y.^2).^2;
        a(1:end - 1) = a(1:end - 1) + 1200 * x(1:end - 1).^2 - 400 * x(2:end) + 2;
        a(2:end) = a(2:end) + 200;
        b = -400 * x(1:end - 1);
        A = spdiags(a, 0, d, d) + sparse(1:d - 1, 2:d, b, d, d) + sparse(2:d, 1:d - 1, b, d, d);
    end

    function v = D(x, p)
        r = x(2:end) - x(1:end - 1).^2;
        dr = p(2:end) - 2 * x(1:end - 1) .* p(1:end - 1) - p(1:end - 1).^2;
        t = atan(x - 1);
        dt = atan2(p, 1 + (x - 1) .* (x + p - 1));
        v = sum(100 * dr .* (2 * r + dr) ...
                + p(1:end - 1) .* (2 * (x(1:end - 1) - 1) + p(1:end - 1))) ...
            + sum(s .* dt .* (2 * t + dt));
    end

end
