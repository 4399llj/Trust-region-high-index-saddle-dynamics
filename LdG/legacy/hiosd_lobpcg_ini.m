function [V, alpha, iter] = hiosd_lobpcg_ini(x,options)
%% reading parameters
n=length(x);
% innp can compute multi - one innp as a column vector
if isfield(options,'innp')
    innp = @(x,y)options.innp(x,y);
    nor = @(x)sqrt(innp(x,x));
else
    disp('HiOSD: standard inner product.');
    innp = @(x,y)(x'*y);
    nor = @(x)sqrt(x'*x);
end

if isfield(options,'H')
    H = @(xx,vv,ll)options.H(xx,vv,ll);
else H = @ (xx,vv,ll) - ( F (xx + ll * vv) - F (xx - ll * vv) ) / (2 * ll);
end

if isfield(options,'k')
%     if isfield(options,'V')
%         if options.k ~= size(options.V,2)
%             options.k = size(options.V,2);
%             disp('HiOSD: k is reset as the column - number of V.');
%         end
%     end
elseif isfield(options,'V')
    options.k = size(options.V,2);
else
    options.k = 1;
    disp('HiOSD: k is reset as 1.');
end
k = options.k;

V = zeros(n,k);
if isfield(options,'V')
    if rank(options.V) < k
        V = randn(n,k);
    else
        V = options.V;
    end
else
    V = randn(n,k);
end

for i = 1 : k
    V(:,i) = V(:,i) - V(:,1:i-1) * innp(V(:,1:i-1), V(:,i));
    V(:,i) = V(:,i)/nor(V(:,i));
end

if isfield(options,'dt')
    dt=options.dt;
else
    dt=0.01;
end

if isfield(options,'l')
    l=options.l;
else
    l=1e-6;
end
if isfield(options,'minl')
    minl=options.minl;
else
    minl=1e-6;
end
if isfield(options,'grad')
     F = options.grad;
end
if isfield(options,'maxiter')
    maxiter=options.maxiter;
else
    maxiter=1e3;
end
if isfield(options,'precond')
    T=options.precond;
else
    T=@(xx)xx;
end
%% setting environment
alpha (1 : k, 1) = 0;
H = @ (xx,vv,ll)  ( F (xx + ll * vv) - F (xx - ll * vv) ) / (2 * ll);
U = zeros(n,k);
W = zeros(n,k);
Y = zeros(n,k);
%% bb iterations
iter = 1;
Vp = zeros(n,k);
Up = zeros(n,k);
while iter<maxiter
    emp = true( 3 * k, 1 );
    for i = 1 : k
        U(:,i) = H (x, V(:,i), l);
        alpha(i) = innp (V(:,i), U(:,i));
        W(:,i) =  ( U(:,i) - alpha(i) * V(:,i) ) ;
        W(:,i) = W(:,i) - V * innp( V, W(:,i) ) ;
        W(:,i) = W(:,i) - W(:,1:i-1) * innp( W(:,1:i-1), W(:,i) ) ;
        nrmW = nor(W(:,i)) ;
        if nrmW > 1e-8
            W(:,i) = W(:,i) / nrmW ;
            Y(:,i) = H (x, W(:,i), l);
        else
            W(:,i) = 0;
            Y(:,i) = 0;
            emp(k+i) = false;
        end
    end
    for i = 1 : k
        Vp(:,i) = Vp(:,i) - [V W Vp(:,1:i-1)] * innp( [V W Vp(:,1:i-1)], Vp(:,i) );
        nrmVpi = nor(Vp(:,i));
        if nrmVpi > 1e-8
            Vp(:,i) = Vp(:,i) / nrmVpi ;
            Vp(:,i) = Vp(:,i) - [V W Vp(:,1:i-1)] * innp( [V W Vp(:,1:i-1)], Vp(:,i) );
            Vp(:,i) = Vp(:,i) / nor(Vp(:,i));
            Up(:,i) = H (x, Vp(:,i), l);
        else
            Vp(:,i) = 0;
            Up(:,i) = 0;
            emp(2*k+i) = false;
        end
    end
    UU = [V,W,Vp]; UU = UU (:,emp);
    YY = [U,Y,Up]; YY = YY (:,emp);
    Pn = innp(UU,YY);  Pn = (Pn+Pn') / 2;
    [eta, alpha] = eig( Pn, 'vector' );
    
    Vp = V; Up = U;
    V = UU * eta(:,1:k);
    alpha = alpha(1:k);
    if mod(iter, 10)==0
        disp(num2str([iter alpha']));
    end
    l = max(l/(1+dt),minl);
    iter = iter+1;
end

end