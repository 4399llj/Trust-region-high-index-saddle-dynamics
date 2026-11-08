function Hv = hessvecnlc_2d(u,v)
%HESSVECNLC_2D Hessian-vector product for the 2D reduced LdG model.
%
% This function computes Hv = Hessian(u) * v, corresponding exactly to
% ngradnlc_2d(u).
%
% Input:
%   u : 2N x 1 vector, u = [q1(:); q2(:)]
%   v : 2N x 1 vector, v = [v1(:); v2(:)]
%
% Output:
%   Hv : 2N x 1 vector

    global I omega h2 N xi

    if length(u) ~= 2*N
        error('Input u must have length 2*N.');
    end

    if length(v) ~= 2*N
        error('Input v must have length 2*N.');
    end

    q1 = reshape(u(1:N),     I, I);
    q2 = reshape(u(N+1:2*N), I, I);

    v1 = reshape(v(1:N),     I, I);
    v2 = reshape(v(N+1:2*N), I, I);

    %% Linearized elastic part
    % Since the boundary values in ngradnlc_2d are fixed constants,
    % their variations are zero in the Hessian-vector product.

    Lv1 = - ( ...
        [zeros(1,I); v1(1:I-1,:)] ...
      + [v1(2:I,:); zeros(1,I)] ...
      + [zeros(I,1), v1(:,1:I-1)] ...
      + [v1(:,2:I), zeros(I,1)] ...
      - 4*v1 ) / h2;

    Lv2 = - ( ...
        [zeros(1,I); v2(1:I-1,:)] ...
      + [v2(2:I,:); zeros(1,I)] ...
      + [zeros(I,1), v2(:,1:I-1)] ...
      + [v2(:,2:I), zeros(I,1)] ...
      - 4*v2 ) / h2;

    %% Linearized bulk part
    %
    % Gradient:
    %   F1 = -Delta_h q1 + omega*(q1^2+q2^2+xi/2)*q1
    %   F2 = -Delta_h q2 + omega*(q1^2+q2^2+xi/2)*q2
    %
    % Hessian-vector product:
    %   Hv1 = -Delta_h v1
    %       + omega*(3q1^2+q2^2+xi/2)*v1
    %       + omega*(2q1q2)*v2
    %
    %   Hv2 = -Delta_h v2
    %       + omega*(2q1q2)*v1
    %       + omega*(q1^2+3q2^2+xi/2)*v2

    Hv1 = Lv1 ...
        + omega * (3*q1.^2 + q2.^2 + xi/2) .* v1 ...
        + omega * (2*q1.*q2) .* v2;

    Hv2 = Lv2 ...
        + omega * (2*q1.*q2) .* v1 ...
        + omega * (q1.^2 + 3*q2.^2 + xi/2) .* v2;

    Hv = [reshape(Hv1,N,1); reshape(Hv2,N,1)];
end