function [Q, Qinv] = QRPlus(l1, l2)
    n = length(l1);
    A = [l1, l2, eye(n)];
    [Q, ~] = qr(A);
    Qinv = inv(Q);
end

