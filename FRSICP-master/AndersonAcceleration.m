classdef AndersonAcceleration
    properties
        m  % 保存的历史步数
        dim  % 变量的维度
        iter  % 当前迭代次数
        col_idx  % 当前历史矩阵的列索引
        current_u  % 当前的变量值
        current_F  % 当前的函数值差异
        prev_dG  % 前几次迭代的变量更新的差分
        prev_dF  % 前几次迭代的函数值更新的差分
        M  % 正则方程矩阵
        theta  % 正则方程的解
        dF_scale  % 函数值差分的缩放因子
    end
    
    methods
        function obj = AndersonAcceleration(m, dim)
            if nargin < 2
                error('输入参数不足。必须提供 m 和 dim。');
            end
            
            assert(m > 0);
            obj.m = m;
            obj.dim = dim;
            obj.current_u = zeros(dim, 1);
            obj.current_F = zeros(dim, 1);
            obj.prev_dG = zeros(dim, m);
            obj.prev_dF = zeros(dim, m);
            obj.M = zeros(m, m);
            obj.theta = zeros(m, 1);
            obj.dF_scale = zeros(m, 1);
            obj.iter = 0;
            obj.col_idx = 1;
        end
        
        function obj = replace(obj, u)
            obj.current_u = u;
        end
        
        function u = compute(obj, g)
            assert(obj.iter >= 0);

            G = g(:);
            obj.current_F = G - obj.current_u;

            if obj.iter == 0
                obj.prev_dF(:, 1) = -obj.current_F;
                obj.prev_dG(:, 1) = -G;
                obj.current_u = G;
            else
                obj.prev_dF(:, obj.col_idx) = obj.prev_dF(:, obj.col_idx) + obj.current_F;
                obj.prev_dG(:, obj.col_idx) = obj.prev_dG(:, obj.col_idx) + G;

                eps = 1e-14;
                scale = max(eps, norm(obj.prev_dF(:, obj.col_idx)));
                obj.dF_scale(obj.col_idx) = scale;
                obj.prev_dF(:, obj.col_idx) = obj.prev_dF(:, obj.col_idx) / scale;

                m_k = min(obj.m, obj.iter);

                if m_k == 1
                    obj.theta(1) = 0;
                    dF_sqrnorm = obj.prev_dF(:, obj.col_idx)' * obj.prev_dF(:, obj.col_idx);
                    obj.M(1, 1) = dF_sqrnorm;
                    dF_norm = sqrt(dF_sqrnorm);

                    if dF_norm > eps
                        obj.theta(1) = (obj.prev_dF(:, obj.col_idx) / dF_norm)' * (obj.current_F / dF_norm);
                    end
                else
                    new_inner_prod = (obj.prev_dF(:, obj.col_idx)' * obj.prev_dF(:, 1:m_k))';
                    obj.M(obj.col_idx, 1:m_k) = new_inner_prod';
                    obj.M(1:m_k, obj.col_idx) = new_inner_prod;

                    % Solve normal equation
                    obj.theta(1:m_k) = obj.M(1:m_k, 1:m_k) \ (obj.prev_dF(:, 1:m_k)' * obj.current_F);
                end

                % Use rescaled theta to compute new u
                obj.current_u = G - obj.prev_dG(:, 1:m_k) * (obj.theta(1:m_k) ./ obj.dF_scale(1:m_k));
                obj.col_idx = mod(obj.col_idx, obj.m) + 1;
                obj.prev_dF(:, obj.col_idx) = -obj.current_F;
                obj.prev_dG(:, obj.col_idx) = -G;
            end

            obj.iter = obj.iter + 1;
            u = obj.current_u;
        end

        function obj = reset(obj, u)
            obj.iter = 0;
            obj.col_idx = 1;
            obj.current_u = u;
        end
    end
end
