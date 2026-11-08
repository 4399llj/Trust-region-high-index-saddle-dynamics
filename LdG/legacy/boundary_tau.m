function tau = boundary_tau(p,d,Delta,innp)
%BOUNDARY_TAU Find tau >= 0 such that ||p + tau d|| = Delta.

    a = innp(d,d);
    b = 2*innp(p,d);
    c = innp(p,p) - Delta^2;

    disc = b^2 - 4*a*c;
    disc = max(disc,0);

    tau = (-b + sqrt(disc))/(2*a);
end