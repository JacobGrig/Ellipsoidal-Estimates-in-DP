function  [lNew, xPlus, XPlus, xMinus, XMinus] = SolveSet(A, B, X1, x1, P, p, t1, t, step, r, l)
    disp(t)  
    n = length(A(t1));
%     expAt = @(split) exp(trapz(split, A(split)));
%     expA_cur = @(tau) exp(trapz(tau : new_step : t1, reshape(A(tau : new_step : t1)));
    %expA_t1_t = expA(t1 - t);
    MFundFunc = @(tau, t) MFund(tau, A, n, t);
    %expA_t_t1 = MFundFunc(t1, t);
    if t ~= t1
        [~, xPlus] = ode45(@(t, x) ForODE(t, x, A, B, p), [t1, t], x1);
    else
        xPlus = x1;
        xMinus = x1;
        XPlus = repmat(X1, 1, 1, r);
        XMinus = repmat(X1, 1, 1, r);
        lNew = MFundFunc(t1, t).'*l(:, 1);
        return
    end
    xPlus = xPlus(end, :);
    xMinus = xPlus(end,:);
    
    splitting = t : step : t1;
    while length(splitting) == 1
        step = step/10;
        splitting = t:step:t1;
    end
%     eye_n = eye(n);
%     l = eye_n(:, 1 : r);
    XPlus = zeros(n, n, r);
    XMinus = zeros(n, n, r);
    for j = 1 : r
        l1 = l(:, j);
        %lt = expA_t1_t.'*l1;
%         p1 = sqrt(dot(l1, X1*l1));
%         p2 = @(tau) sqrt(dot(l1,MFundFunc(t1, tau)*B(tau)*P(tau)*(B(tau).')*(MFundFunc(t1, tau).')*l1));
%         p2Vect = zeros(1, length(splitting));
%         intMatr = zeros(n, n, length(splitting));
%         for i = 1:length(splitting)
%             k = splitting(i);
%             p2Vect(i) = p2(k);
%             intMatr(:, :, i) = MFundFunc(t, k)*B(k)*P(k)*(B(k).')*(MFundFunc(t, k).') ./ p2(k);
%         end
%         XPlus(:, :, j) = (p1 + trapz(splitting, p2Vect))*(expA_t_t1 * X1 * (expA_t_t1.') / p1...
%             + trapz(splitting, intMatr, 3));
        XPlus(:, :, j) = Plus(t, A, X1, P, B, n, l1, t1);
        XMinus(:, :, j) = Minus(t, A, X1, P, B, n, l1, t1);
    end
    lNew = MFundFunc(t1, t).'*l(:, 1);
end

function res = ForPlus(t, x, A, X1, P, B, n, l1, t1)
    matr = reshape(x, n, n);
    MF = MFund(t1, A, n, t);
    Pi = sqrt(dot(l1,MF*B(t)*P(t)*(B(t).')*(MF.')*l1))...
        /sqrt(dot(l1, MF*matr*(MF.')*l1));
    res = -Pi*matr - Pi^(-1)*B(t)*P(t)*(B(t).') + A(t)*matr + matr*(A(t).');
    res = res(:);
end

function res = Plus(t, A, X1, P, B, n, l1, t1)
    res = X1;
    start = res(:);
    if (t ~= t1)
        [~, traect] = ode45(@(t,x) ForPlus(t, x, A, X1, P, B, n, l1, t1), [t1, t], start);
        res = reshape(traect(end, :), n, n);
    end
end

function res = ForOOO(t, x, A, X1, P, B, n, l1, t1)
    matr = reshape(x, n, n);
    a = sqrtm(P(t))*(B(t).')*(MFund(t1,A, n, t).')*l1;
    am = sqrtm(P(t))*(B(t).')*(MFund(t1,A, n, t).');
%     b = sqrt(dot(l1,MFund(t1,A, n, t)*B(t)*P(t)*B(t).'*MFund(t1,A, n, t).'*l1))...
%        /sqrt(dot(l1, X1*l1))*sqrtm(X1)*l1;
    b = norm(a)/norm(sqrtm(X1)*l1)*sqrtm(X1);
    S = b/am;
    %Qb = B(t)*P(t)*(B(t).');
    %H = sqrtm(matr)\S*Qb;
    res = A(t) * matr + matr * (A(t).') - S*sqrtm(P(t))*(B(t).');
    %res = A(t) * matr + matr * (A(t).') - (H.') * matr - matr * H;
    res = (res.')*res;
    res = res(:);
end

function res = Minus(t, A, X1, P, B, n, l1, t1)
    res = sqrtm(X1);
    start = res(:);
    if (t ~= t1)
        [~, traect] = ode45(@(t,x) ForOOO(t, x, A, X1, P, B, n, l1, t1), [t1, t], start);
        res = reshape(traect(end, :), n, n);
    end
end

function res = ForODE(t, x, A, B, p)
    res = A(t)*x + B(t)*p(t);
end

function res = ode_fund(t, fund, A_matr, n)
    fund_matr = reshape(fund, n, n);
    res_matr = A_matr(t) * fund_matr;
    res = res_matr(:);
end

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