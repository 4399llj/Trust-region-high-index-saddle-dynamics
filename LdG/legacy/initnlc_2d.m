clear;
addpath ./nlc_2d
addpath ./nlcdata100
global I N omega x y h h2 xi E0 s
h = 1/50;  x = -1 : h : 1;   y = -1 : h : 1;   I = length(x) - 2;   N = I ^ 2;
%parameters
B = 0.64 * 10^4;  C = 0.35 * 10^4;  
xi = -B^2/2/C/C; E0 = xi^2/8 * 0 ;s = sqrt(2*abs(xi));h2=h^2;
G=2;inpf=@(x,y)x'*(G.*y)*h2;
omega =  50;
E=@(uu)enenlc_2d(uu); F=@(uu)ngradnlc_2d(uu);   

ldnm = @(NM)load(sprintf('./nlcdata%d/L_%.0f_%s',I+1, omega, NM));
dp = @(NM,op)drawpath_beta(NM,op)+drawpath_2d(NM,op);

optionsini = hiosdoptions(1);
optionsini.dt = 0;
optionsini.innpfunc=inpf;
optionsini.outputX = 0; 
optionsini.outputp = 100;
optionsini.outputd = 0;
optionsini.maxiter = 500;

options = hiosdoptions(1);
options.dt = 0.0002;
options.betat = 5;       options.betau = 0.2;
options.gammamax = 0.1;  options.gammamin = 0.001;
options.innpfunc = inpf; options.epsf = 1e-6;
options.outputX = 50;    options.maxiter = 1000000;
options.outputp = 1000;
options.outputd = 500;   options.draws='drawnlc_2d(x);';
return;

options.dt = 0.00002;
options.gammamax = 10;  options.gammamin = 0.1;

NAME='N+1';
u=ldnm(NAME);u=u(1:N*2);drawnlc_2d(u);drawnow;
optionsini.k=6;rng(0);optionsini.V=randn(length(u),optionsini.k);
[V, alpha] = hiosd_lobpcg(F,u,optionsini);

V(:,[4 3])=V(:,[3 4])*[-V(N+7401,3) V(N+7401,4);-V(N+7401,4) -V(N+7401,3)]/norm(V(N+7401,[3 4]));

options.k=2;options.V=V(:,[1 2]);
[u1, output1] = hiosd_sirqit(F,u+0.00001*V(:,2),options);

options.k=0;options.V=V(:,[]);
[u1, output1] = hiosd_sirqit(F,u+0.00001*V(:,1),options);
[u0, output0] = hiosd_sirqit(F,u-0.00001*V(:,1),options);

options.k=1;options.V=V(:,1)-V(:,2);
[u1, output1] = hiosd_sirqit(F,u+0.00001*(V(:,1)+V(:,2)),options);


options.k=3;options.V=[V(:,[1 2]) V(:,3)-V(:,4)];
[u1, output1] = hiosd_sirqit(F,u+0.00001*(V(:,3)+V(:,4)),options);