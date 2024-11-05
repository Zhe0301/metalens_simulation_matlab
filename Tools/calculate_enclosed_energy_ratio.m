function [radius, energy_ratio] = calculate_enclosed_energy_ratio(I, G, center_x, center_y,varargin)
% 根据质心位置计算圈入能量比
% I: 光强矩阵
% G: 坐标网路类
% center_x, center_y: 圈入的中心位置
% 可选参数
% max_radius: 圈入最大半径，默认为网格长度的一半
% total_energy: 输入光的总能量，默认梯形积分计算
% min_points: 模拟圆的最小细分点数，默认为128
    p = inputParser;    
    addParameter(p,'max_radius',0); 
    addParameter(p,'total_energy',0);
    addParameter(p,'min_points',128);
    parse(p,varargin{:});       
    max_radius = p.Results.max_radius;
    total_energy = p.Results.total_energy;
    min_points = p.Results.min_points;
    % 设置参数默认值
    if max_radius == 0
        max_radius = G.length/2;
    end
    if total_energy == 0
        total_energy = trapz(G.axis,trapz(G.axis,I,2));
    end
        
    % 计算x和y坐标的最大范围
    x_max_dist = max(abs(G.axis - center_x));
    y_max_dist = max(abs(G.axis - center_y));
    radius_limit = min(x_max_dist, y_max_dist);  % 最大有效半径
    if max_radius>radius_limit
        fprintf("invalid max_radius, radius_limit=%.4f mm",radius_limit)
        max_radius = radius_limit;
    end
        

    % 筛选需要计算的区域
    x_idx = G.axis >= (center_x - radius_limit) & G.axis <= (center_x + max_radius);
    y_idx = G.axis >= (center_y - radius_limit) & G.axis <= (center_y + max_radius);
    x_c = G.axis(x_idx);
    y_c = G.axis(y_idx);
    y_c = y_c';
    I_c = I(y_idx,x_idx);

    % 创建插值对象，用于计算非整数位置的光强
    interpI = griddedInterpolant({x_c, y_c}, I_c, 'linear', 'none');

    % 初始化半径和能量比数组
    radius = 0:G.step:max_radius;
    energy_ratio = zeros(size(radius));

    % 计算每个半径下的圈入能量比
    for i = 1:length(radius)
        if i == 1
            energy_ratio(i) = 0;
        end
        radii = radius(i);
        
        if i > min_points
            interp_points = i;
        else
            interp_points = min_points;  
        end

        x_interp = linspace(center_x - radii, center_x + radii, interp_points);
        y_interp = linspace(center_y - radii, center_y + radii, interp_points);
        I_c_interp = interpI({x_interp,y_interp'});
        [X_interp,Y_interp] = meshgrid(x_interp,y_interp);
        distances = sqrt((X_interp - center_x).^2 + (Y_interp - center_y).^2);
        I_c_interp(distances>radii) = 0;
        encircled_energy = trapz(y_interp,trapz(x_interp,I_c_interp,2));
        energy_ratio(i) = encircled_energy / total_energy;
    end

end
