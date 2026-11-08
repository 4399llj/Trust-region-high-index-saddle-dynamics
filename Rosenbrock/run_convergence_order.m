function summary = run_convergence_order()
    % Section 4.1 convergence-order test: d=100, index 3, sneg=-5000.
    here = fileparts(mfilename('fullpath'));
    out = fullfile(here, 'results', 'order');
    if ~exist(out, 'dir')
        mkdir(out);
    end

    %% parameter setting
    m = rosenbrock_model(100, 3, -5000);
    x0 = 1.1 * ones(100, 1);
    rows = {};

    %% convergence-order test
    for method = {'CP', 'N'}
        o = struct('method', method{1}, 'dt', 1e-5, 'exactSpectrum', true, ...
                 'tol', 1e-5, 'maxiter', 1e5, 'delta', 1e-3);
        rng(20261107);
        r = paper_solve(m, x0, o);
        actual_index = sum(eig(full(m.hessian(r.x))) < 0);
        assert(r.converged && norm(r.x - m.exact) < 1e-4 && actual_index == 3);

        %% estimate local convergence order
        error = r.error(:);
        error_n = error(1:end - 1);
        error_next = error(2:end);
        local = error_n < 1e-2 & error_next > 1e-13 & error_next < error_n;
        local_samples = sum(local);
        estimated_order = NaN;

        % Use the same local fitting window as plot_convergence_order.py.
        if local_samples >= 3
            coefficient = polyfit(log10(error_n(local)), ...
                                  log10(error_next(local)), 1);
            estimated_order = coefficient(1);
        else
            warning('Too few local samples to estimate the order for %s.', method{1});
        end

        if strcmp(method{1}, 'CP')
            reference_order = 1;
        else
            reference_order = 2;
        end

        order = struct('estimated_order', estimated_order, ...
                       'reference_order', reference_order, ...
                       'local_samples', local_samples, 'local_mask', local);
        save(fullfile(out, ['d100_' method{1} '.mat']), 'r', 'x0', 'order');
        rows(end + 1, :) = {method{1}, r.iterations, r.residual(end), ...
                           error(end), actual_index, estimated_order, ...
                           reference_order, local_samples};

        error_pairs = table((0:numel(error_n) - 1)', error_n, error_next, ...
                            error_next ./ error_n, error_next ./ error_n.^2, local, ...
                            'VariableNames', {'n', 'error_n', 'error_next', ...
                            'linear_ratio', 'quadratic_ratio', 'local_fit_sample'});
        writetable(error_pairs, fullfile(out, [method{1} '_error_pairs.csv']));
        fprintf('TRHiSD-%s: estimated order = %.6f, reference = %d, local samples = %d\n', ...
                method{1}, estimated_order, reference_order, local_samples);
    end

    %% output data
    summary = cell2table(rows, 'VariableNames', ...
                        {'method', 'iterations', 'residual', 'error', 'index', ...
                         'estimated_order', 'reference_order', 'local_samples'});
    writetable(summary, fullfile(out, 'summary.csv'));
    disp(summary);
end
