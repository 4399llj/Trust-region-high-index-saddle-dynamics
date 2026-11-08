function run_reproduction(dimensions)
    % Run all four methods with explicit common initial states.
    % Default: both dimensions. run_reproduction(2) is a quick low-dimensional run.
    if nargin == 0
        dimensions = [2, 1000];
    end
    assert(all(ismember(dimensions, [2, 1000])));
    here = fileparts(mfilename('fullpath'));
    out = fullfile(here, 'results');
    if ~exist(out, 'dir')
        mkdir(out);
    end

    %% parameter setting
    template = load(fullfile(here, 'data', 'plot_data.mat'));
    D2 = template.D2;
    D1000 = template.D1000;
    methods = {'CP', 'BB', 'HiSD', 'AHiSD'};
    rows = {};

    %% simulation
    for d = dimensions
        if d == 2
            k = 1;
            s = -100;
            beta = 1e-3;
        else
            k = 5;
            s = -50000;
            beta = 1e-5;
        end
        m = rosenbrock_model(d, k, s);
        x0 = 1.1 * ones(d, 1);
        runs = cell(1, 4);
        for j = 1:4
            o = struct('method', methods{j}, 'dt', beta, 'bbmin', beta / 5, 'bbmax', 10 * beta, ...
                     'tol', 1e-5, 'maxiter', 1e5, 'eigtol', 1e-3, 'delta', 1e-3, 'momentum', 0.3);
            rng(20261107);
            r = paper_solve(m, x0, o);
            runs{j} = r;
            assert(numel(r.error) == r.iterations + 1 && all(isfinite(r.error)));
            actual_index = sum(eig(full(m.hessian(r.x))) < 0);
            rows(end + 1, :) = {d, methods{j}, r.iterations, r.converged, ...
                               r.residual(end), r.error(end), actual_index, r.seconds};
            if d == 2
                D2.lines(j).X = 0:r.iterations;
                D2.lines(j).Y = r.error(:)';
            else
                D1000.lines(j).X = 0:r.iterations;
                D1000.lines(j).Y = r.error(:)';
            end
            fprintf('d=%d %s: iterations=%d converged=%d error=%.12g\n', ...
                    d, methods{j}, r.iterations, r.converged, r.error(end));
        end
        save(fullfile(out, sprintf('d%d_all_methods.mat', d)), 'runs', 'x0', 'methods');
    end

    %% output data
    summary = cell2table(rows, 'VariableNames', ...
                        {'dimension', 'method', 'iterations', 'converged', ...
                         'residual', 'error', 'index', 'seconds'});
    writetable(summary, fullfile(out, 'summary.csv'));
    % Partial runs are for checking the solver and do not create a mixed-data plot.
    if all(ismember([2, 1000], dimensions))
        save(fullfile(out, 'plot_data.mat'), 'D2', 'D1000');
    end
end
