function res = MFund(tau, AMatr,n, t)
    eye_matr = eye(n);
    eye_lin = eye_matr(:);
    if (tau == t)
        res = eye(n);
    else
        [~, res_vect] = ode45(@(t, fund) ode_fund(t, fund, AMatr, n), [tau, t], eye_lin);
        res = reshape(res_vect(end, :), n, n);
    end
end