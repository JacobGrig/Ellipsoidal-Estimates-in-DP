function border_points = DrawEllips(XPlus, xPlus, l1, l2, epsilon1, epsilon2)
    n = length(xPlus);
    l1 = l1/norm(l1, 2);
    l2 = l2 - dot(l1, l2)*l1;
    l2 = l2/norm(l2,2);
    xPlusProj = l1*dot(xPlus, l1) + l2*dot(xPlus, l2);
    proj_dir = xPlus - xPlusProj;
    if (norm(proj_dir) < epsilon1)
        [Q, ~] = QRPlus(l1, l2);
        l3 = Q(:, 3);
        xPlus = xPlus + epsilon2 * l3;
        xPlusProj = l1*dot(xPlus, l1) + l2*dot(xPlus, l2);
        proj_dir = xPlus - xPlusProj;
    end
    proj_dir = proj_dir / norm(proj_dir);
    splitting = 0 : epsilon2 : 2*pi;
    border_points = zeros(n, length(splitting));
    for i = 1 : length(splitting)
        alpha = splitting(i);
        gamma = 1.0;
        beta = 1.0;
        init_point = xPlus;
        while (beta >= epsilon1) || (abs(gamma) >= epsilon1)
            gamma = 1.0;
            beta = 1.0;
            while true
                new_point = init_point + beta * (l1 * cos(alpha) + l2 * sin(alpha));
                if isOnBorderOfEllips(new_point, xPlus, XPlus, epsilon2)
                    break
                else
                    if isInEllips(new_point, xPlus, XPlus, epsilon2)
                        init_point = new_point;
                    else
                        beta = beta / 2;
                    end
                end
            end
            % new point is on the border in direction alpha

            new_point_pos = new_point + proj_dir;
            new_point_neg = new_point - proj_dir;

            gamma = isCloserToEllips(new_point_pos, new_point_neg, xPlus, XPlus);
            
            saved_point = new_point;
            init_point = new_point;
            while true%abs(gamma) >= epsilon
                new_point = init_point + gamma * proj_dir;
                if isOnBorderOfEllips(new_point, xPlus, XPlus, epsilon2)
                    break
                else
                    if isInEllips(new_point, xPlus, XPlus, epsilon2)
                        init_point = new_point;
                    else
                        gamma = gamma / 2;
                    end
                end
            end
            if abs(gamma) < epsilon1
                %break
            else
                init_point = (saved_point + new_point) / 2;
            end
        end
        border_points(:, i) = l1*dot(init_point, l1) + l2*dot(init_point, l2);
    end
end

function res = isInEllips(point, xPlus, XPlus, epsilon)
    res = (dot(point - xPlus, XPlus * (point - xPlus)) <= 1 + epsilon);
end

function res = isOnBorderOfEllips(point, xPlus, XPlus, epsilon)
    res = ((abs(dot(point - xPlus, XPlus * (point - xPlus)) - 1)) <= epsilon);
end

function res = isCloserToEllips(point_pos, point_neg, xPlus, XPlus)
    res = 2*((abs(dot(point_pos - xPlus, XPlus * (point_pos - xPlus)) - 1)) <= ...
             (abs(dot(point_neg - xPlus, XPlus * (point_neg - xPlus)) - 1))) - 1;
end