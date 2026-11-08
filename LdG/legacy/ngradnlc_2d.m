function F = ngradnlc_2d(u)
global I omega h2 s N xi
q1 = reshape( u(    1:  N), I, I );
q2 = reshape( u(  N+1:2*N), I, I );

otrQ2xi = omega * ( (q1.^2 + q2.^2) + xi / 2);

LB1 = - ( [ones(1,I)*(-s/2);q1(1:I-1,:)] + [q1(2:I,:);ones(1,I)*(-s/2)] + [ones(I,1)*(s/2) q1(:,1:I-1)] + [q1(:,2:I) ones(I,1)*(s/2)] - 4 * q1) /h2  ...
    + otrQ2xi .* q1 ;
LB2 = - ( [zeros(1,I) ;q2(1:I-1,:)] + [q2(2:I,:);zeros(1,I) ] + [zeros(I,1)  q2(:,1:I-1)] + [q2(:,2:I) zeros(I,1) ] - 4 * q2) /h2  ...
    + otrQ2xi .* q2 ;
F = [reshape(LB1,N,1);reshape(LB2,N,1)];
end