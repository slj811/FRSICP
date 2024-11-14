function T = align_symm(p1, p2, n1, n2, w)
    w = w / sum(w);  % 归一化权重
    w2 = sqrt(w);
    meanp1 = mean(p1, 2);
    meanp2 = mean(p2, 2);
    p1 = p1 - meanp1;
    p2 = p2 - meanp2;
    n = n1 + n2;

    % 初始化矩阵和向量
    A = zeros(6, 6);
    b = zeros(6, 1);

    for i = 1:size(p1, 2)
        c = cross(p1(:, i), n(:, i));
        d = p2(:, i) - p1(:, i);
        x = [c; n(:, i)];
        A = A + (x * x') * w(i);
        b = b + x * dot(d, n(:, i)) * w(i);
       
    end

    x = A \ b;  % 解线性方程组

    rot = x(1:3);
    trans = x(4:6);
    rotangle = norm(rot);
    TR = rotation_matrix(rotangle, rot);
    trans = trans + meanp2 - TR(1:3, 1:3) * meanp1;
    T = [TR(1:3, 1:3), trans; 0, 0, 0, 1];
end
