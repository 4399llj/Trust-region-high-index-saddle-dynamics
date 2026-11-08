function E = enenlc_2d(u)
global I omega E0 h2 N s xi
q1 = reshape( u(    1:  N), I, I );
q2 = reshape( u(  N+1:2*N), I, I );

Ee=[norm([ s/2*ones(I,1) q1]-[q1  s/2*ones(I,1)],'fro')^2;
    norm([    zeros(I,1) q2]-[q2     zeros(I,1)],'fro')^2;
    norm([-s/2*ones(1,I);q1]-[q1;-s/2*ones(1,I)],'fro')^2;
    norm([    zeros(1,I);q2]-[q2;    zeros(1,I)],'fro')^2;];

trQ2 = 2*(q1.^2 + q2.^2);
% Eb = (trQ2 + xi).^2 / 8;
Eb = trQ2.^2/8 + xi/4*trQ2;% + E0;

E = sum(Ee) + omega * sum(sum(Eb)) * h2;
end