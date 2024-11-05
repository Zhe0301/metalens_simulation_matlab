function [x_interp,y_interp,I_interp] = extract_spot(I, G,center_x, center_y,max_radius,interp_factor)
% 提取光斑，并适当插值
% 输入参数
% I
    
    % 计算x和y坐标的最大范围
    x_max_dist = max(abs(G.axis - center_x));
    y_max_dist = max(abs(G.axis - center_y));
    radius_limit = min(x_max_dist, y_max_dist);  % 最大有效半径
    if max_radius>radius_limit
        fprintf("invalid max_radius, radius_limit=%.4f mm",radius_limit)
        max_radius = radius_limit;
    end
        

    % 筛选需要计算的区域
    x_idx = G.axis >= (center_x - max_radius) & G.axis <= (center_x + max_radius);
    y_idx = G.axis >= (center_y - max_radius) & G.axis <= (center_y + max_radius);
    x_c = G.axis(x_idx);
    y_c = G.axis(y_idx);
    y_c = y_c';
    I_c = I(y_idx,x_idx);
     % 创建插值对象，用于计算非整数位置的光强
    interpI = griddedInterpolant({x_c, y_c}, I_c, 'linear', 'none');
    x_interp = linspace(max(x_c), min(y_c), length(y_c)*interp_factor);
    y_interp = linspace(max(x_c), min(y_c), length(y_c)*interp_factor);
    y_interp = y_interp';
    I_interp = interpI({x_interp,y_interp});
end