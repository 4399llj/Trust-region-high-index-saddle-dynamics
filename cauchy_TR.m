function [s,H_f,nor_f] = cauchy_TR(f,phi,delta,flag,H,tor,innp)
%CAUCHY_TR  cauchy point method for trust region
% flag control the problem to solve min or max
% flag == 1 f descend direction => s descend step size > 0 
% flag == 2 f ascend direction => s ascend step size > 0  
    nor_f = sqrt(innp(f,f));
    if nor_f > tor/2
        H_f = innp(f,H(phi , f));
        if H_f*flag > 0
            s = min([delta/nor_f,flag*nor_f^2/H_f]);
        else
            s = delta/nor_f;
        end
    else
        H_f = 0;
        nor_f = 0;
        s = 0;
    end
end
