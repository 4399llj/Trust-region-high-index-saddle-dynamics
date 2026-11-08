function [phi, options, meta] = ldg_setup(state, k)
    % Self-contained initialization for the two-component 2D LdG model.
    if nargin < 1
        state = 'T1';
    end
    if nargin < 2
        k = 6;
    end

    %% model parameters
    global I N omega x y h h2 xi E0 s
    h = 1 / 50;
    x = -1:h:1;
    y = -1:h:1;
    I = length(x) - 2;
    N = I^2;
    h2 = h^2;
    B = .64e4;
    C = .35e4;
    xi = -B^2 / (2 * C * C);
    E0 = 0;
    s = sqrt(2 * abs(xi));
    omega = 50;

    %% initial state
    assert(ischar(state) || isstring(state));
    state = char(state);
    assert(~isempty(regexp(state, '^[A-Za-z0-9+_-]+$', 'once')), 'Invalid state name');
    here = fileparts(mfilename('fullpath'));
    source = fullfile(here, 'data', 'nlcdata100', ['L_50_' state]);
    phi = load(source, '-ascii');
    assert(isvector(phi) && numel(phi) >= 2 * N);
    phi = phi(1:2 * N);
    phi = phi(:);
    assert(k >= 1 && k == floor(k) && k < 2 * N);

    %% solver options
    options = osdoptions(1);
    options.output = 100;
    options.maxiter = 100000;
    options.dt = 1e-4;
    options.tor = 1e-5;
    options.k = k;
    options.innp = @(a, b) a' * (2 * b) * h2;
    options.grad = @ngradnlc_2d;
    options.E = @enenlc_2d;
    options.dif_E = @dif_enenlc_2d;
    % Solver calls H(u,v); eigen-initializer calls H(u,v,l).
    options.H = @(u, v, varargin) hessvecnlc_2d(u, v);
    rng(20261107);
    options.V = randn(2 * N, k);
    meta = struct('state', state, 'k', k, 'h', h, 'I', I, ...
                  'omega', omega, 'xi', xi, 'seed', 20261107, 'source', source);
end
