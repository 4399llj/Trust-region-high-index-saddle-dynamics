clear;
addpath ./nlc
addpath ./nlcdata100
global I N s lambda x y h A B C E0 A2C A4C B2C h2
h = 1/50;  x = -1 : h : 1;   y = -1 : h : 1;   I = length(x) - 2;   N = I ^ 2;
%parameters
B = 0.64 * 10^4;  C = 0.35 * 10^4;   A = -B^2/3/C;   
s = ( B + sqrt(B^2 + 24 * abs(A) * C) )/4/C;
E0=s^2*A/(6*C)-s^3*B/(27*C)+s^4/18;
A2C=A/2/C; A4C=A/4/C; B2C=B/2/C; h2=h^2;
G=[2*ones(N,1);2*ones(N,1);6*ones(N,1);2*ones(N,1);2*ones(N,1)];
inpf=@(x,y)x'*(G.*y)*h2;

lambda =     15  ;

E=@(uu)enenlc(uu);
F=@(uu)ngradnlc(uu);   
H = @ (xx,vv,ll) - ( F (xx + ll * vv) - F (xx - ll * vv) ) / (2 * ll);

options=osdoptions(1);options.outputX=100; options.dt=0.0002;
options.gammamax=0.1;options.gammamin=0.001;options.betat=5;options.betau=0.2;
options.innpfunc=inpf;   options.epsf=1e-6;
options.outputp = 100; 
options.draws = "drawnlc(x)"; options.outputd = 100;


ldnm = @(NM)load(sprintf('./nlcdata%d/L_%d_%s',I+1, lambda, NM));
savepath=@(NUM)(['drawpath(' num2str(NUM) ',output' num2str(NUM) ');output' num2str(NUM) '.X=[];']);

%%
NAME='WORS';
u=ldnm(NAME);drawnlc(u);drawnow;
options.k=3;options.maxiter=500;options.V=randn(5*N,options.k);
[V, alpha] = hiosd_lobpcg_ini(F,u,options);

options.k=1;options.V=V(:,1);options.outputX=10;options.maxiter=100000;
[u1, output1] = hiosd_sirqit(F,u+0.001*V(:,2),options);
options.k=1;options.V=V(:,1);options.outputX=10;options.maxiter=100000;
[u2, output2] = hiosd_sirqit(F,u-0.001*V(:,2),options);

NAME='BDX';
u=ldnm(NAME);drawnlc(u);drawnow;
options.k=2;options.maxiter=500;options.V=randn(5*N,options.k);
[V, alpha] = hiosd_lobpcg_ini(F,u,options);

options.k=0;options.V=V(:,[]);options.outputX=10;options.maxiter=100000;
[u3, output3] = hiosd_sirqit(F,u+0.001*V(:,1),options);
[u4, output4] = hiosd_sirqit(F,u-0.001*V(:,1),options);

NAME='BDY';
u=ldnm(NAME);drawnlc(u);drawnow;
options.k=2;options.maxiter=500;options.V=randn(5*N,options.k);
[V, alpha] = hiosd_lobpcg_ini(F,u,options);

options.k=0;options.V=V(:,[]);options.outputX=10;options.maxiter=100000;
[u5, output5] = hiosd_sirqit(F,u+0.001*V(:,1),options);
[u6, output6] = hiosd_sirqit(F,u-0.001*V(:,1),options);

for NUM=9:20
    eval(savepath(NUM));
end 