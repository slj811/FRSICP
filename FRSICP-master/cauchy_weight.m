function w = cauchy_weight(r, u)
    % 计算权重
    w = 1 ./ (1 + (r / u).^2);
end