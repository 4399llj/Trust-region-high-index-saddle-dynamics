clear;
addpath ./nlc
addpath ./nlcdata100
mkdir video
global I N s lambda x y h A B C E0 A2C A4C B2C h2
h = 1/50;  x = -1 : h : 1;   y = -1 : h : 1;   I = length(x) - 2;   N = I ^ 2;
%parameters
B = 0.64 * 10^4;  C = 0.35 * 10^4;   A = -B^2/3/C;   
s = ( B + sqrt(B^2 + 24 * abs(A) * C) )/4/C;
E0=s^2*A/(6*C)-s^3*B/(27*C)+s^4/18;
A2C=A/2/C; A4C=A/4/C; B2C=B/2/C; h2=h^2;
G=[2*ones(N,1);2*ones(N,1);6*ones(N,1);2*ones(N,1);2*ones(N,1)];
inpf=@(x,y)x'*(G.*y)*h2;

lambda =     50   ;

E=@(uu)enenlc(uu);
F=@(uu)ngradnlc(uu);   
H = @ (xx,vv,ll) - ( F (xx + ll * vv) - F (xx - ll * vv) ) / (2 * ll);

options=osdoptions(1);options.outputX=100; options.dt=0.0002;
options.gammamax=0.1;options.gammamin=0.02;options.betat=5;options.betau=0.2;
options.innpfunc=inpf;   options.epsf=1e-5;
options.outputp = 100; 
options.draws = "drawnlc(x)"; options.outputd = 100;

ldnm = @(NM)load(sprintf('./nlcdata%d/L_%d_%s',I+1, lambda, NM));
%% WORS=>S1
NAME='WORS';u=ldnm(NAME);drawnlc(u,1,NAME);
options.k=7;    options.maxiter=500;   options.V=randn(5*N,options.k);
[V, alpha] = hiosd_lobpcg_ini(F,u,options);
V(:,6)=V(:,6)*sign(V(N+N/2+0.5,6));

options.k=4;    options.V=V(:,1:4);    options.maxiter=100000;
[u1, output1] = hiosd_lobpcg(F,u+0.0001*V(:,6),options);
drawpath(1,output1);output1.X=[];

%% S1=>T+1
NAME='S1';u=ldnm(NAME);drawnlc(u,1,NAME);
options.k=5;    options.maxiter=500;   options.V=randn(5*N,options.k);
[V, alpha] = hiosd_lobpcg_ini(F,u,options);
V(:,[3 4])=[-V(:,[3 4])*V(16301,[3 4])'/norm(V(16301,[3 4])) -V(:,[3 4])*V(13233,[3 4])'/norm(V(13233,[3 4]))];
V(:,5)=sign(V(16564,5))*V(:,5);

options.k=3;    options.V=V(:,[1 2 4]) ;    options.maxiter=100000;
[u2, output2] = hiosd_lobpcg(F,u+0.0001*V(:,3),options);
drawpath(2,output2);output2.X=[];
%% S1=>T-2
options.k=3;    options.V=V(:,1:3);    options.maxiter=100000;
[u3, output3] = hiosd_lobpcg(F,u+0.0001*V(:,4),options);
drawpath(3,output3);output3.X=[];
%% S1=>S+1
options.k=5;    options.V=V(:,1:5);    options.maxiter=100000;
[u4, output4] = hiosd_lobpcg(F,u+0.0001*V(:,5),options);
drawpath(4,output4);output4.X=[];
%% S1=>S-1
options.k=5;    options.V=V(:,1:5);    options.maxiter=100000;
[u5, output5] = hiosd_lobpcg(F,u-0.0001*V(:,5),options);
drawpath(5,output5);output5.X=[];

%% S+1=>C+
NAME='S+1';u=ldnm(NAME);drawnlc(u,1,NAME);
options.k=5;    options.maxiter=500;   options.V=randn(5*N,options.k);
[V, alpha] = hiosd_lobpcg_ini(F,u,options);
V(:,5)=-V(:,5)*sign(V(N+N/2+0.5,5));

options.k=2;    options.V=V(:,1:2);    options.maxiter=100000;
[u6, output6] = hiosd_lobpcg(F,u+0.0001*V(:,5),options);
drawpath(6,output6);output6.X=[];

%% C+=>M+E
NAME='C+';u=ldnm(NAME);drawnlc(u,1,NAME);
options.k=4;    options.maxiter=500;   options.V=randn(5*N,options.k);
[V, alpha] = hiosd_lobpcg_ini(F,u,options);
V(:,3)=-V(:,[3 4])*V(17276,[3 4])'/norm(V(17276,[3 4]));

options.k=3;    options.V=V(:,1:3);    options.maxiter=100000;
[u7, output7] = hiosd_lobpcg(F,u+0.0001*V(:,3),options);
drawpath(7,output7);output7.X=[];

%% M+E=>N+3
NAME='M+E';u=ldnm(NAME);drawnlc(u,1,NAME);
options.k=5;    options.maxiter=500;   options.V=randn(5*N,options.k);
[V, alpha] = hiosd_lobpcg_ini(F,u,options);
V(:,4)=-sign(V(14715,4))*V(:,4);

options.k=4;    options.V=V(:,1:4);    options.maxiter=5000000; options.outputX=1000;options.epsf=1e-10;
[u8, output8] = hiosd_lobpcg(F,u+0.0001*V(:,4),options);
drawpath(8,output8);output8.X=[];


