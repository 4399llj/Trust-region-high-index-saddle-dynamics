function v_eps = pre_eps(x0,x1,V)
    % give two states and predict the best eps of V
    % x0 initial state
    % x1 final state
    % V the eigenvector of x0
    v_eps=zeros(size(V,2),1);
    V1=V*v_eps+x0-x1;
    H=V'*V;
    f=H\V'*V1;
    while norm(f)>1e-4
        v_eps=v_eps-f;
        V1=V*v_eps+x0-x1;
        f=H\V'*V1;
        disp(norm(V*v_eps+x0-x1))
    end
end