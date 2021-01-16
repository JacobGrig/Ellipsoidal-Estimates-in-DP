function [new_l1, new_l2, new_border_points] = drawProj(border_points, l1, l2)
    [Q, Qinv] = QRPlus(l1, l2);
    new_l1 = Q(:, 1);
    new_l2 = Q(:, 2);
    new_border_points = Qinv * border_points;
end