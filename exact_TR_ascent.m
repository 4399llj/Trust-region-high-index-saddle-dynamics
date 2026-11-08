function [p, pred, info] = exact_TR_ascent(f1, V, alpha, Delta, innp, nor)
%EXACT_TR_ASCENT Exact trust-region Newton ascent step on the unstable subspace.
%
% This function solves
%
%     max_p  f1^T p + 1/2 p^T H p
%     s.t.   p in span(V),  ||p|| <= Delta,
%
% where V contains the negative eigenspace basis and alpha contains the
% corresponding Hessian eigenvalues.
%
% If H v_i = alpha_i v_i with alpha_i < 0, then on span(V),
%
%     A = -H|_V
%
% is positive definite with eigenvalues a_i = -alpha_i. The subproblem becomes
%
%     max_z  fV^T z - 1/2 z^T A z
%     s.t.   ||V z|| <= Delta.
%
% Input:
%   f1    : projected gradient on V, f1 = V*(V^T f)
%   V     : N x k basis of the unstable subspace
%   alpha : 1 x k or k x 1 Hessian eigenvalues on V, normally alpha_i < 0
%   Delta : trust-region radius
%   innp  : inner product function
%   nor   : norm function induced by innp
%
% Output:
%   p     : exact trust-region ascent step in span(V)
%   pred  : predicted increase m(p)-m(0)
%   info  : diagnostic information



    info = struct();

    N = size(V,1);
    k = size(V,2);

    if Delta <= 0
        p = zeros(N,1);
        pred = 0;
        info.sigma = Inf;
        info.boundary = true;
        info.z = zeros(k,1);
        return;
    end

    alpha = alpha(:);

    if length(alpha) ~= k
        error('Length of alpha must equal the number of columns of V.');
    end

    % Coefficients of f1 in the basis V.
    fV = innp(V, f1);
    fV = real(fV(:));

    if norm(fV) == 0
        p = zeros(N,1);
        pred = 0;
        info.sigma = 0;
        info.boundary = false;
        info.z = zeros(k,1);
        return;
    end

    % A = -H|_V should be positive definite.
    a = -real(alpha);


    % Unconstrained Newton ascent step:
    %     z_N = A^{-1} fV.
    zN = fV ./ a;
    pN = V * zN;

    if nor(pN) <= Delta
        z = zN;
        p = pN;
        sigma = 0;
        boundary = false;
    else
        % Boundary case:
        % Find sigma > 0 such that
        %     || V * (fV ./ (a + sigma)) || = Delta.
        sigma_low = 0;
        sigma_high = 1;

        z_high = fV ./ (a + sigma_high);
        p_high = V * z_high;

        while nor(p_high) > Delta
            sigma_high = 2 * sigma_high;
            z_high = fV ./ (a + sigma_high);
            p_high = V * z_high;
        end

        for it = 1:80
            sigma_mid = 0.5 * (sigma_low + sigma_high);

            z_mid = fV ./ (a + sigma_mid);
            p_mid = V * z_mid;

            if nor(p_mid) > Delta
                sigma_low = sigma_mid;
            else
                sigma_high = sigma_mid;
            end
        end

        sigma = sigma_high;
        z = fV ./ (a + sigma);
        p = V * z;
        boundary = true;
    end

    % Predicted increase:
    %
    %     pred = fV^T z + 1/2 z^T diag(alpha) z
    %          = fV^T z - 1/2 z^T diag(a) z.
    %
    pred = sum(fV .* z) + 0.5 * sum(alpha .* (z.^2));
    pred = real(pred);

    if abs(pred) <= 1e-14
        pred = 0;
    end

    info.sigma = sigma;
    info.boundary = boundary;
    info.z = z;
    info.normp = nor(p);
    info.pred = pred;
end