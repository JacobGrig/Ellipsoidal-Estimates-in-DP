function res = PreSet(l, Q, q)
    res = q + l/(sqrt(dot(l, Q*l)));
end

