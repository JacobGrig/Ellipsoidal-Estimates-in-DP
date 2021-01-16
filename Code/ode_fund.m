function res = ode_fund(t, fund, A_matr, n)
    fund_matr = reshape(fund, n, n);
    res_matr = A_matr(t) * fund_matr;
    res = res_matr(:);
end