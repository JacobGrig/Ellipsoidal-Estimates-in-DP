%% Input
clc;
n = 3;
m = 3;
r = 10;
A = @(x) [cos(x) -sin(x) 0;
          sin(x) cos(x) 0;
          0 0 1];
B = @(x)[0 cos(x) 0;
         0 1 1;
         1 0 1];
t1 = 1;
t0 = 0.2;
X1 = [1 0 0;
      0 1 0;
      0 0 1];
x1 = [0; 3; 2];
P = @(x) [sin(x) 0 0;
          0 1 0;
          0 0 1];
p = @(x) [0; 0; 0];
step = 0.01;
step1 = 0.2;

%% Basic Input
clc;
n = 3;
m = 3;
r = 60;
A = @(x) [1 0 0;
          0 1 0;
          0 0 1];
B = @(x)[1 0 0;
         0 0 0;
         0 0 0] + eye(3)*step;
t1 = 1;
t0 = 0;
X1 = [1 0 0;
      0 1 0;
      0 0 1];
x1 = [0; 0; 0];
P = @(x) [1 0 0;
          0 1 0;
          0 0 1];
p = @(x) [0; 0; 0];
step = 0.01;
step1 = 0.2;

%% Execution
clc;

l = rand(n, r);

% l = [1 0 0;
%      0 1 0;
%      0 0 1;
%      1 0 1;
%      1 1 1];
% l = l.';
splitting = t1: -step1 : t0;
xPlus = zeros(n, length(splitting));
XPlus = zeros(n, n, r, length(splitting));
xMinus = zeros(n, length(splitting));
XMinus = zeros(n, n, r, length(splitting));
lNew = zeros(n, length(splitting));
for i = 1:length(splitting)
    [lNew(:, i), xPlus(:, i), XPlus(:, :, :, i),xMinus(:, i), XMinus(:, :, :, i)] ...
            = SolveSet(A, B, X1, x1, P, p, t1, splitting(i), step, r, l);    
end

%% Draw Dynamic
l1 = lNew(:, end);
l2 = l(:, 2);
epsilon1 = 1e-5;
epsilon2 = 1e-2;
len = length(0 : epsilon2 : 2*pi);
BorderPointsMatrix = zeros(2, len, length(splitting), r);
BorderPointsMatrixMinus = zeros(2, len, length(splitting), r);

for j = 1:r
    disp(j)
    for i = 1 : length(splitting)
        %l1 = lNew(:, i);
        border_points = DrawEllips(XPlus(:, :, j, i), xPlus(:,i), l1, l2, epsilon1, epsilon2);
        disp(i)
        [new_l1, new_l2, new_border_points] = drawProj(border_points, l1, l2);
        BorderPointsMatrix(:, :, i, j) = new_border_points(1:2, :);
    end
end

for j = 1:r
    disp(j)
    for i = 1 : length(splitting)
        %l1 = lNew(:, i);
        border_points = DrawEllips(XMinus(:, :, j, i), xMinus(:,i), l1, l2, epsilon1, epsilon2);
        disp(i)
        [new_l1, new_l2, new_border_points] = drawProj(border_points, l1, l2);
        BorderPointsMatrixMinus(:, :, i, j) = new_border_points(1:2, :);
    end
end

%% Plots
clc;
MatrRes = cell(1,length(splitting));
for j = 1:length(splitting)
    xa = BorderPointsMatrix(1, :, j, 1);
    ya = BorderPointsMatrix(2, :, j, 1);
    for i = 2:r
        [xa, ya] = polybool('intersection', xa, ya, ...
        BorderPointsMatrix(1, :, j, i), BorderPointsMatrix(2, :, j, i));
    end
    MatrRes(j) = {[xa; ya]};
end

MatrResMinus = cell(1,length(splitting));
for j = 1:length(splitting)
    xa = BorderPointsMatrixMinus(1, :, j, 1);
    ya = BorderPointsMatrixMinus(2, :, j, 1);
    for i = 2:r
        [xa, ya] = polybool('union', xa, ya, ...
        BorderPointsMatrix(1, :, j, i), BorderPointsMatrixMinus(2, :, j, i));
    end
    MatrResMinus(j) = {[xa; ya]};
end

%%
clc;
max = 0;
for i = 1:length(MatrResMinus)
    Matr = MatrResMinus{i};
    if length(Matr(1, :)) > max
        max = length(Matr(1, :));
    end
end


XMin = zeros(max, length(splitting));
YMinus = zeros(max, length(splitting));
ZMinus = zeros(max, length(splitting));

for i = 1:length(MatrResMinus)
    Matr = MatrResMinus{i};
    if length(Matr(1, :)) < max
        XMin(1:length(Matr(1, :)), i) = Matr(1, :).';
        YMinus(1:length(Matr(1, :)), i) = Matr(2, :).';
        XMin(length(Matr(1, :)) + 1 : end, i) = ones(length(length(Matr(1, :)) + 1 : max), 1)*Matr(1, end);
        YMinus(length(Matr(1, :)) + 1 : end, i) = ones(length(length(Matr(1, :)) + 1 : max), 1)*Matr(2, end);
    else
        XMin(:, i) = Matr(1,:).';
        YMinus(:, i) = Matr(2,:).';
    end
    
    ZMinus(:,i) = ones(max, 1)*splitting(i);
end

max = 0;
for i = 1:length(MatrRes)
    Matr = MatrRes{i};
    if length(Matr(1, :)) > max
        max = length(Matr(1, :));
    end
end


X = zeros(max, length(splitting));
Y = zeros(max, length(splitting));
Z = zeros(max, length(splitting));

for i = 1:length(MatrRes)
    Matr = MatrRes{i};
    if length(Matr(1, :)) < max
        X(1:length(Matr(1, :)), i) = Matr(1, :).';
        Y(1:length(Matr(1, :)), i) = Matr(2, :).';
        X(length(Matr(1, :)) + 1 : end, i) = ones(length(length(Matr(1, :)) + 1 : max), 1)*Matr(1, end);
        Y(length(Matr(1, :)) + 1 : end, i) = ones(length(length(Matr(1, :)) + 1 : max), 1)*Matr(2, end);
    else
        X(:, i) = Matr(1,:).';
        Y(:, i) = Matr(2,:).';
    end
    
    Z(:,i) = ones(max, 1)*splitting(i);
end

%%
surf(Z, X, Y,'FaceAlpha',0.5, 'EdgeColor', 'none');
hold on;
s = surf(ZMinus, XMin, YMinus,'FaceAlpha',0.5, 'EdgeColor', 'none');