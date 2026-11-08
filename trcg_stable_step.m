function p = trcg_stable_step(phi, g2, V, Delta, H, innp, nor)
%TRCG_STABLE_STEP Steihaug--Toint truncated CG on the stable subspace W.
%
% It approximately solves
%
%   min_p  g2' p + 1/2 p' H p
%   s.t.   p in W, ||p|| <= Delta,
%
% where W is the orthogonal complement of span(V).
%
% Input:
%   phi   : current base point
%   g2    : projected gradient Q g
%   V     : basis of the unstable subspace, assumed orthonormal
%   Delta : trust-region radius on W
%   H     : Hessian-vector product, H(phi,v)
%   innp  : inner product
%   nor   : norm induced by innp
%   cgopt : struct with fields tol and maxit
%
% Output:
%   p     : actual descent step in W, to be added to phi

    if Delta <= 0
        p = zeros(size(g2));
        return;
    end

    tol = 1e-4;
    maxit = 75;


    % Projection onto W
    projW = @(z) z - V*innp(V,z);

    g2 = projW(g2);
    ng = nor(g2);

    if ng == 0
        p = zeros(size(g2));
        return;
    end

    % Initial point p = 0
    p = zeros(size(g2));

    % Residual of model gradient: r = g2 + H p.
    % At p=0, r=g2.
    r = g2;

    % CG direction
    d = -r;
    d = projW(d);

    rr = innp(r,r);
    rr0 = rr;

    for it = 1:maxit
        Hd = H(phi,d);
        Hd = projW(Hd);

        dHd = innp(d,Hd);

        % Negative curvature: go to the boundary
        if dHd <= 0
            tau = boundary_tau(p,d,Delta,innp);
            p = p + tau*d;
            p = projW(p);
            return;
        end

        alpha_cg = rr / dHd;

        % Boundary hit
        if nor(p + alpha_cg*d) >= Delta
            tau = boundary_tau(p,d,Delta,innp);
            p = p + tau*d;
            p = projW(p);
            return;
        end

        % Standard CG step
        p = p + alpha_cg*d;
        r_new = r + alpha_cg*Hd;
        r_new = projW(r_new);

        rr_new = innp(r_new,r_new);

        % Stopping condition
        if sqrt(rr_new) <= tol * sqrt(rr0)
            p = projW(p);
            return;
        end

        beta_cg = rr_new / rr;
        d = -r_new + beta_cg*d;
        d = projW(d);

        r = r_new;
        rr = rr_new;
    end

    p = projW(p);
end