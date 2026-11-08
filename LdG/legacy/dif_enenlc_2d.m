function dE = dif_enenlc_2d(u,p)
%DIF_ENENLC_2D Compute E(u+p)-E(u) for enenlc_2d without subtracting
% two large energy values.
%
% This function matches the energy code enenlc_2d(u):
%
%   E = sum(Ee) + omega * sum(sum(Eb)) * h2;
%
% where
%   trQ2 = 2*(q1.^2 + q2.^2);
%   Eb   = trQ2.^2/8 + xi/4*trQ2.
%
% Input:
%   u : 2N x 1 vector, u = [q1(:); q2(:)]
%   p : 2N x 1 increment vector
%
% Output:
%   dE = E(u+p)-E(u)

    global I omega h2 N s xi

    if length(u) ~= 2*N
        error('Input u must have length 2*N.');
    end
    if length(p) ~= 2*N
        error('Input p must have length 2*N.');
    end

    q1  = reshape(u(1:N),     I, I);
    q2  = reshape(u(N+1:2*N), I, I);

    dq1 = reshape(p(1:N),     I, I);
    dq2 = reshape(p(N+1:2*N), I, I);

    %% --------------------------------------------------------------------
    % 1. Elastic energy difference
    %
    % For each difference term A, we compute
    %
    %   ||A+dA||_F^2 - ||A||_F^2
    %   = 2*sum(A.*dA) + sum(dA.^2).
    %
    % Boundary values are fixed constants, so their increments are zero.
    %% --------------------------------------------------------------------

    % q1 horizontal differences
    A1  = [ s/2*ones(I,1), q1] - [q1, s/2*ones(I,1)];
    dA1 = [ zeros(I,1), dq1]   - [dq1, zeros(I,1)];

    % q2 horizontal differences
    A2  = [ zeros(I,1), q2] - [q2, zeros(I,1)];
    dA2 = [ zeros(I,1), dq2] - [dq2, zeros(I,1)];

    % q1 vertical differences
    A3  = [-s/2*ones(1,I); q1] - [q1; -s/2*ones(1,I)];
    dA3 = [ zeros(1,I); dq1]   - [dq1; zeros(1,I)];

    % q2 vertical differences
    A4  = [ zeros(1,I); q2] - [q2; zeros(1,I)];
    dA4 = [ zeros(1,I); dq2] - [dq2; zeros(1,I)];

    dEe1 = 2*sum(A1(:).*dA1(:)) + sum(dA1(:).^2);
    dEe2 = 2*sum(A2(:).*dA2(:)) + sum(dA2(:).^2);
    dEe3 = 2*sum(A3(:).*dA3(:)) + sum(dA3(:).^2);
    dEe4 = 2*sum(A4(:).*dA4(:)) + sum(dA4(:).^2);

    dE_elastic = dEe1 + dEe2 + dEe3 + dEe4;

    %% --------------------------------------------------------------------
    % 2. Bulk energy difference
    %
    % trQ2 = 2*(q1.^2 + q2.^2)
    %
    % Eb = trQ2.^2/8 + xi/4*trQ2.
    %
    % Let dtrQ2 = trQ2(u+p)-trQ2(u). Then
    %
    % Eb(u+p)-Eb(u)
    % =
    %   1/8 * dtrQ2 .* (2*trQ2 + dtrQ2)
    %   + xi/4 * dtrQ2.
    %% --------------------------------------------------------------------

    trQ2 = 2*(q1.^2 + q2.^2);

    dtrQ2 = ...
        4*(q1.*dq1 + q2.*dq2) ...
        + 2*(dq1.^2 + dq2.^2);

    dEb = ...
        (1/8) * dtrQ2 .* (2*trQ2 + dtrQ2) ...
        + (xi/4) * dtrQ2;

    dE_bulk = omega * sum(dEb(:)) * h2;

    %% --------------------------------------------------------------------
    % 3. Total energy difference
    %% --------------------------------------------------------------------

    dE = dE_elastic + dE_bulk;
end