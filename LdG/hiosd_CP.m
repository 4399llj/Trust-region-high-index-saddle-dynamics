function [phi,V,alpha,iter] = hiosd_CP(phi,options)
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
    end
    if isfield(options,'output')
        output = options.output;
    end
    if isfield(options,'dif_E')
        dif_E = options.dif_E;
    end
    if isfield(options,"k")
        k = options.k;
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
    if isfield(options,"H")
        H = options.H;
    else
        ll = 1e-8;
        H = @(uu,vv) -(grad(uu+ll*vv)-grad(uu-ll*vv))/(2*ll);
    end
    if isfield(options,"maxiter")
        maxiter = options.maxiter;
    end
    if isfield(options,"tor")
        tor = options.tor;
    else
        tor = 1e-5;
    end
    if isfield(options,"proj")
        Proj = options.proj;
    end
   %%  environment setting
    iter = 0;
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
    delta = [0.001,0.001]; % trust region ascend & descend
    gamma_1 = -1;
    gamma_2 = -1;

%% iteration process
while iter <= maxiter
    %% LOBPCG
    emp = true( 3 * k, 1 );
    for i = 1 : k
        V(:,i) = V(:,i);
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
    Vp = V;
    V = UU * real(eta(:,1:k)); 
    alpha = real(alpha(1:k))';
    %% trust region - ascend direction -exact calculation
    f = grad(phi);
    f_1 = V*innp(V,f);
    f_2 = f - f_1;
    while 1
        %% trust region - cauchy point
        if gamma_1 <= 0
            % % cauchy point method
            [s_1,H_f1,norf_1] = cauchy_TR(f_1,phi,delta(1),-1,H,tor,innp);
            p_1 = s_1*f_1;            
        end
        if  gamma_2 <= 0
            [s_2,H_f2,norf_2] = cauchy_TR(f_2,phi,delta(2),1,H,tor,innp);
            p_2 = s_2*f_2;
        end

        %% update trust region domain size

        gamma_1 = (dif_E(phi,p_1))/(s_1*norf_1^2+s_1^2*H_f1/2)';

        if gamma_1 == inf
            gamma_1 = 1;
        elseif gamma_1 == -inf
            gamma_1 = -1;
        elseif isnan(gamma_1)
            gamma_1 = 0.5;
        end

        gamma_2 = (dif_E(phi,p_2))/(s_2*norf_2^2+s_2^2*H_f2/2)';

        if gamma_2 == inf
            gamma_2 = 1;
        elseif gamma_2 == -inf
            gamma_2 = -1;
        elseif isnan(gamma_2)
            gamma_2 = 0.5;
        end

        delta(1) = update_domain(gamma_1,delta(1));
        delta(2) = update_domain(gamma_2,delta(2));

        if gamma_1 > 0 && gamma_2 > 0
            phi = phi + p_1 - p_2;
            gamma_1 = -1;
            gamma_2 = -1;
            break
        end
    end
    %% output data
    if mod(iter,output) == 0
       
       disp([nor(f)])
       % disp([t0 alpha E(phi_A) ])
       % disp([nor(grad(phi_A)) ])
    end    
    norf = nor(f);
    if norf < tor
        break
    end
    
    iter = iter + 1;
end

