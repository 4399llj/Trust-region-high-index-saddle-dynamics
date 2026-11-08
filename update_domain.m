function delta = update_domain(gamma,delta)
% UPDATE_DOMAIN 此处显示有关此函数的摘要
%   此处显示详细说明
    % common update
    if gamma> 0.75
        delta = 2*delta;
    elseif gamma <= 0.25
        delta = delta/4;
    end
    % stable setting
    if delta >= 1e-2
        delta = 1e-2;
    end

end

