classdef Grid < handle
    % 建立二维正方形仿真网格 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % num_points (int): 网格在一个维度上的点数
    % length (int): 网格在一个维度上的长度，单位mm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Last modified: 2024.10.29
    % Author: zwz
    properties
        num_points
        length
        step
        step_fft
        half_num_points
        vector
        axis
        axis_fft
        d2_x
        d2_y
        d2_square
        d2_r
        d2_theta
        d2_fft_x
        d2_fft_y
    end
    
    methods
        function G = Grid(num_points,length)
            G.num_points = num_points;
            G.length = length;
            G.step = G.length / G.num_points; % 离散实空间步长
            G.step_fft = 1 / G.length; % 离散频率步长
            G.half_num_points = G.num_points / 2;
            
            G.vector = 1:G.num_points;
            G.axis = -G.length / 2 + G.step / 2 + (G.vector - 1) * G.step; % 离散后实空间坐标轴
            G.axis_fft = -1 / (2 * G.step) + (G.vector - 1) / G.length; % 离散后频谱空间坐标轴
            
            [G.d2_x, G.d2_y] = meshgrid(G.axis, G.axis); % 离散后的坐标矩阵
            G.d2_square = G.d2_x .^ 2 + G.d2_y .^ 2;
            G.d2_r = sqrt(G.d2_square); % 离散后的极坐标半径
            G.d2_theta = atan2(G.d2_y, G.d2_x); % 离散后的极坐标极角，atan2可自动判断象限，值域为-pi~pi；atan值域为-pi/2~pi/2
            
            [G.d2_fft_x, G.d2_fft_y] = meshgrid(G.axis_fft, G.axis_fft);
        end
    end
end
