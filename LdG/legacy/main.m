global I N omega x y h h2 xi E0 s
h = 1/50;  x = -1 : h : 1;   y = -1 : h : 1;   I = length(x) - 2;   N = I ^ 2;
%parameters
B = 0.64 * 10^4;  C = 0.35 * 10^4;  
xi = -B^2/2/C/C; E0 = xi^2/8 * 0 ;s = sqrt(2*abs(xi));h2=h^2;
ll = 1e-15;
G=2;
omega =  50;


ldnm = @(NM)load(sprintf('./nlcdata%d/L_%.0f_%s',I+1, omega, NM));
dp = @(NM,op)drawpath_beta(NM,op)+drawpath_2d(NM,op);

%% polymers model energy and F=-grad E;Energy
options = osdoptions(1); options.output  = 100;
grad  = @(uu) ngradnlc_2d(uu);
% H = @ (xx,vv)  ( grad(xx + ll * vv) - grad(xx - ll * vv) ) / (2 * ll);
H = @ (xx,vv) hessvecnlc_2d(xx,vv);
innp  = @(x,y)x'*(G.*y)*h2;
E = @(uu)enenlc_2d(uu);
dif_E = @(uu,pp) dif_enenlc_2d(uu,pp);
options.innp = innp;
options.grad = grad;
options.E = E;
options.H = H;
options.maxiter = 100000;
options.dt = 1e-4;
options.tor = 1e-5; % torelate error
options.dif_E = dif_E;
%% from u0_multipod
NAME='T1';
phi=ldnm(NAME);

phi = phi(1:N*2);
options.k=6; 
options.V=rand(N*2,options.k);
% options.V = V_L_50_N1(:,1:3);
% [phi_3,V_3,alpha_3,store_3] = hiosd(L_50_BD_+0.1*V_L_50_BD_(:,4),options);
% [L_50_M1,V_L_50_M1,alpha_1,iter] = hiosd_CP(L_50_M1,options);
%[phi_3,~,~,iter] = hiosd_CPCG(L_50_N1+0.1*V_L_50_N1(:,4),options);
[V_L_50_T, ~, iter] = hiosd_lobpcg_ini(phi,options);
