function [q1, q2] = vec2q_2d(u)
global I s N
q1 = zeros(I+2,I+2);
q2 = zeros(I+2,I+2);

q1(1  , 2:I+1) = -s/2;
q1(I+2, 2:I+1) = -s/2;
q1(2:I+1, 1  ) =  s/2;
q1(2:I+1, I+2) =  s/2;

q1([1 I+2],[1 I+2])=0;

q1(2:I+1, 2:I+1) = reshape( u(    1:  N), I, I );
q2(2:I+1, 2:I+1) = reshape( u(  N+1:2*N), I, I );
end