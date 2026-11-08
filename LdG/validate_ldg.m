function checks = validate_ldg()
    % Model and dependency check; this does not claim a full landscape search.
    [u, o] = ldg_setup('T1', 6);
    rng(20261107);
    v = randn(size(u));
    v = v / norm(v);
    eps = 1e-5;
    g = o.grad(u);
    hv = o.H(u, v);

    %% derivative checks
    checks.gradient = abs((o.E(u + eps * v) - o.E(u - eps * v)) / (2 * eps) ...
                          - o.innp(g, v)) / max(1, abs(o.innp(g, v)));
    checks.hessian = norm((o.grad(u + eps * v) - o.grad(u - eps * v)) / (2 * eps) ...
                         - hv) / max(1, norm(hv));
    p = eps * v;
    checks.energy_difference = abs(o.dif_E(u, p) - (o.E(u + p) - o.E(u))) ...
                               / max(1, abs(o.dif_E(u, p)));
    assert(all(isfinite(g)) && all(isfinite(hv)));
    assert(checks.gradient < 1e-5 && checks.hessian < 1e-5 && checks.energy_difference < 1e-8);
    assert(norm(o.H(u, v, 1e-6) - hv) == 0);
    % Exercise the formerly incompatible 3-argument Hessian callback.
    o.maxiter = 2;
    [V, eigenvalues, ~] = hiosd_lobpcg_ini(u, o);
    assert(all(isfinite(V(:))) && all(isfinite(eigenvalues)));
    checks.orthogonality = norm(o.innp(V, V) - eye(6), 'fro');
    assert(checks.orthogonality < 1e-6);
    here = fileparts(mfilename('fullpath'));
    out = fullfile(here, 'results');
    if ~exist(out, 'dir')
        mkdir(out);
    end
    save(fullfile(out, 'validation.mat'), 'checks');
    disp(checks);
end
