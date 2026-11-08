function result = run_ldg(method, state, k, maxiter)
    % Example: run_ldg('CP','T1',6,100000).
    % Match state/index/perturbation to the intended landscape transition.
    if nargin < 1
        method = 'eigenspace';
    end
    if nargin < 2
        state = 'T1';
    end
    if nargin < 3
        k = 6;
    end
    if nargin < 4
        maxiter = 100000;
    end

    %% parameter setting
    [phi, options, meta] = ldg_setup(state, k);
    options.maxiter = maxiter;
    initial_phi = phi;
    t = tic;

    %% iteration process
    switch upper(method)
        case 'EIGENSPACE'
            [V, alpha, iter] = hiosd_lobpcg_ini(phi, options);
        case 'CP'
            [phi, V, alpha, iter] = hiosd_CP(phi, options);
        case {'CPCG', 'TCG'}
            [phi, V, alpha, iter] = hiosd_CPCG(phi, options);
        case 'HISD'
            [phi, V, alpha, iter] = hiosd(phi, options);
        otherwise
            error('Unknown method: %s', method);
    end

    %% output data
    g = options.grad(phi);
    residual = sqrt(options.innp(g, g));
    result = struct('phi', phi, 'initial_phi', initial_phi, ...
                    'V', V, 'alpha', alpha, 'iter', iter, ...
                  'seconds', toc(t), 'energy', options.E(phi), 'residual', residual, ...
                    'converged', ~strcmpi(method, 'eigenspace') && residual < options.tor, ...
                    'meta', meta, 'method', method);
    here = fileparts(mfilename('fullpath'));
    out = fullfile(here, 'results');
    if ~exist(out, 'dir')
        mkdir(out);
    end
    save(fullfile(out, sprintf('%s_%s_k%d.mat', upper(method), state, k)), 'result', 'options');
    fprintf('%s: iter=%d energy=%.12g residual=%.12g converged=%d\n', ...
            method, iter, result.energy, residual, result.converged);
end
