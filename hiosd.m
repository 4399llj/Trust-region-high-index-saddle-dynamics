function [phi,V,alpha,iter] = hiosd(phi,options)
    %% parameter setting
    N = length(phi);
    if isfield(options,'innp')
        innp = options.innp;
        nor = @(x)sqrt(innp(x,x));
    else
        disp('HiOSD: standard inner product.');
        innp = @(x,y)(x'*y);
        nor = @(x)sqrt(innp(x,x));
    end
    if isfield(options,'grad')
        grad = options.grad;
    else
        grad = @(uu)grad_AC(uu);
    end
    if isfield(options,'output')
        output = options.output;
    end
    if isfield(options,'dt')
        dt = options.dt;
    else
        dt = 1e-4;
    end
    if isfield(options,'E')
        E = options.E;
    else
        E = @(uu)energy(uu);
    end
    if isfield(options,"k")
        if isfield(options,"V")
            k = options.k;
        end
    else
        disp('HiOSD: k is reset to 1');
        k = 1;
    end
    if isfield(options,'V')
        if rank(options.V) < k
            V = randn(N,k);
        else
            V = options.V;
        end
    else
        V = randn(N,k);
    end
    % data is not save by default
    if isfield(options,'save')
        save = options.save;
    else
        save = -1;
    end
    if isfield(options,"H")
        H = options.H;
    else
        H = @(uu,vv,ll) -(grad(uu+ll*vv)-grad(uu-ll*vv))/(2*ll);
    end
    if isfield(options,"maxiter")
        maxiter = options.maxiter;
    end
    if isfield(options,"tor")
        tor = options.tor;
    else
        tor = 1e-5;
    end
   %%  environment setting
    iter = 0;
    store =[];
    alpha = zeros(1,k);
    for i = 1:k
        V(:,i) = V(:,i) - V(:,1:i-1) * innp(V(:,i),V(:,1:i-1))';
        V(:,i) = V(:,i)/nor(V(:,i));
    end
    % H = @ (xx,vv,ll) - ( grad(xx + ll * vv) - grad(xx - ll * vv) ) / (2 * ll);
    Vp = zeros(N,k);
    Up = zeros(N,k);
    U = zeros(N,k);
    W = zeros(N,k);
    Y = zeros(N,k);
    f = grad(phi);
    norf = nor(f);
    gp = zeros(N,1);
    betat = 10; % bb parameter
    betau = 0.2; % bb parameter
    Dphi = f - 2 * V * innp(V,f); Dphi = Dphi * dt;
    adp_e = 0;
    %% iteration process
while iter <= maxiter
    if norf < tor
        break
    end
    emp = true( 3 * k, 1 );
    
    g = f - 2 * V * innp(V,f);
    Dg = g - gp;
    gp = g;
    beta = abs ( innp(Dphi,Dg) / innp(Dg,Dg) );
    beta = min ( beta, betat * dt);
    beta = max ( beta, betau * dt);
    Dphi = - g * beta;
    phi = phi + Dphi;

    %% LOBPCG
    adp_e = adp_e + norm(Dphi,inf);
    if  adp_e > 1e-3 
        for i = 1 : k
            U(:,i) = H (phi, V(:,i));
            alpha(i) = innp (V(:,i), U(:,i));
            W(:,i) =  ( U(:,i) - alpha(i) * V(:,i) ) ;
            W(:,i) = W(:,i) - V * innp(W(:,i) ,V)' ;
            W(:,i) = W(:,i) - W(:,1:i-1) * innp(W(:,i), W(:,1:i-1))' ;
            nrmW = nor(W(:,i)) ;
            if nrmW > 1e-8
                W(:,i) = W(:,i) / nrmW ;
                Y(:,i) = H (phi, W(:,i));
            else
                % Y(:,i) = H (x, W(:,i), l);
                W(:,i) = 0;
                Y(:,i) = 0;
                emp(k+i) = false;
            end
        end
        for i = 1 : k
            zeta = innp( Vp(:,i),[V W Vp(:,1:i-1)]);
            Vp(:,i) = Vp(:,i) - [V W Vp(:,1:i-1)] * zeta';
            % Up(:,i) = Up(:,i) - [U Y Up(:,1:i-1)] * zeta;
            nrmVpi = nor(Vp(:,i));
            if nrmVpi > 1e-8
                Vp(:,i) = Vp(:,i) / nrmVpi ;
                Vp(:,i) = Vp(:,i) - [V W Vp(:,1:i-1)] * innp(Vp(:,i),[V W Vp(:,1:i-1)])';
                Vp(:,i) = Vp(:,i) / nor(Vp(:,i)) ;
                %Up(:,i) = Up(:,i) / nrmVpi ;
                Up(:,i) = H (phi, Vp(:,i));
            else
                %Up(:,i) = H (x, Vp(:,i), l);
                Vp(:,i) = 0;
                Up(:,i) = 0;
                emp(2*k+i) = false;
            end
        end
        UU = [V,W,Vp]; UU = UU (:,emp); 
        YY = [U,Y,Up]; YY = YY (:,emp);
    
        Pn = innp(UU,YY);  Pn = (Pn+Pn') / 2;
        [eta, alpha] = eig( Pn, 'vector' );
        Vp=V;
        V = UU * real(eta(:,1:k)); 
        alpha = real(alpha(1:k))';
        adp_e = 0;
    end
    %% next grad
    f = grad(phi);
    norf = nor(grad(phi));
    %% output data
    if mod(iter,output) == 0
        disp(norf)
        % disp([t0 alpha E(phi_A) ])
        % disp([nor(grad(phi_A)) ])
    end
    
    [~,sort_index] = sort(alpha);
    V = V(:,sort_index);
    alpha = alpha(sort_index);
    iter = iter+1;
end

