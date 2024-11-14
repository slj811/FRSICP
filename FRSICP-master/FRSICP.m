 function [T, count] = FRSICP(SP,TP,SN,TN, num_components)
[coeff_SP, SP_pca, ~] = pca(SP');
[coeff_TP, TP_pca, ~] = pca(TP');
SP_pca = SP_pca(:, 1:num_components)';
TP_pca = TP_pca(:, 1:num_components)';
    
SN_pca = coeff_SP(:, 1:num_components)' * SN;
TN_pca = coeff_TP(:, 1:num_components)' * TN;

T = eye(4);
Btree = KDTreeSearcher(TP_pca');
[~, dist] = knnsearch(Btree, TP_pca', 'k', 7);
dist = dist(:, 2:7);
u2 = median(median(dist, 2)) / (3 * sqrt(3));
[p1, p2, n1, n2, r] = match_points(SP_pca, TP_pca, SN_pca, TN_pca, Btree);
u1 = median(r);
stop1 = 0; count = 0;
max_iterations = 100;
m = 5;  % Anderson历史步数
dim = numel(T); 
anderson = AndersonAcceleration(m, dim);
%last_energy = inf;
while(~stop1 && count < max_iterations)
   for i = 1:max_iterations - count 
        w = cauchy_weight(r, u1);
        T0 = symm_po_pl(p1, p2, n1, n2, w);
        % Anderson加速
        anderson = anderson.replace(T0(:)); 
        T_accel = reshape(anderson.compute(T0(:)), 4, 4);  
        
        stop2 = norm(T0 - eye(4));
        stop2_accel = norm(T_accel - eye(4));
        
        if stop2_accel < stop2 
            T = T_accel * T;
        else
            T = T0 * T; 
        end
        %p12_T0 = T0 * [SP; ones(1, size(SP, 2))];
        %p12_T_accel = T_accel * [SP; ones(1, size(SP, 2))];
    
        %p1_T0 = p12_T0(1:3, :);
        %p1_T_accel = p12_T_accel(1:3, :);
        %[p1_T0, p2_TP, ~, ~, ~] = match_points(p1_T0, TP, n1, TN, Btree);
        %[p1_T_accel, ~, ~, ~, ~] = match_points(p1_T_accel, TP, n1, TN, Btree);
        
        %rmse_T0 = sqrt(sum(sum((p1_T0 - p2_TP).^2)) / size(p1_T0, 2));
        %rmse_T_accel = sqrt(sum(sum((p1_T_accel - p2_TP).^2)) / size(p1_T_accel, 2));
        
        %if rmse_T_accel < rmse_T0
            %T = T_accel * T;
        %else
            %T = T0 * T; 
        %end
        p12 = T * [SP_pca; ones(1, size(SP_pca, 2))]; 
        p1 = p12(1:3, :);
        n1 = T(1:3, 1:3) * SN_pca;
        [p1, p2, n1, n2, r] = match_points(p1, TP_pca, n1, TN_pca, Btree);
        
        % 检查收敛性
        stop2 = norm(T_accel - eye(4));
        if stop2 < 1e-5
            break;
        end
    end
    if abs(u1 - u2) < 1e-6 || count >= max_iterations
        stop1 = 1;
    end
    count = count + i;
    u1 = u1 / 4;
    if u1 < u2
        u1 = u2;
    end
    
end

p12 = T * [SP_pca; ones(1, size(SP_pca, 2))]; 
p1 = p12(1:3, :); 
n1 = T(1:3, 1:3) * SN_pca;
[idx, dist] = knnsearch(Btree, p1');
inliers = dist < 3 * u2;
p1 = p1(:, inliers); 
n1 = n1(:, inliers); 
p2 = TP_pca(:, idx(inliers)); 
n2 = TN_pca(:, idx(inliers)); 
ww = cauchy_weight(dist(inliers), 3 * u2)';

T0 = symm_po_pl(p1, p2, n1, n2, ww);

T = T0 * T;
T = [coeff_TP, zeros(3,1); zeros(1,3), 1] * T * inv([coeff_SP, zeros(3,1); zeros(1,3), 1]);
mean_SP = mean(SP, 2);
mean_TP = mean(TP, 2);
T(1:3, 4) = T(1:3, 4) + mean_TP - T(1:3, 1:3) * mean_SP;
