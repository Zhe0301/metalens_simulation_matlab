classdef PropOperator < handle
    % Grid网格类的实例，
    % wavelength_vacuum: 真空中波长，
    % dist: 传输距离，
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 可选变量
    % refractive_index: 折射率，默认值为1
    % method: 衍射计算方法，AS为角谱，BL-AS为带限角谱，FFT-DI是基于FFT的直接积分法，默认值为AS
    % paraxial是否进行近轴近似，仅在AS方法下有效
    % gpu_acceleration: 是否使用GPU进行加速计算
    % BL-AS带限角谱理论 ：Band-Limited Angular Spectrum Method for Numerical Simulation of Free-Space Propagation in Far and Near Fields
    % 避免角谱中由于传递函数的采样问题，产生的严重数值误差，考虑倏逝波  
    % FFT-DI 基于FFT的瑞利-索末菲直接积分法：Fast-Fourier-transform based numerical integration method for the Rayleigh–Sommerfeld diffraction formula
    properties
        e_out
        dist
        n
        wavelength
        k_prop
        method
        refractive_index
        mat
        mat_DI
        quality
        gpu_acceleration
    end
    
    methods
        function P = PropOperator(Grid, wavelength_vacuum, dist, varargin)

            p = inputParser;
            addParameter(p,'refractive_index',1); 
            addParameter(p,'method','AS'); 
            addParameter(p,'paraxial',false); 
            addParameter(p,'gpu_acceleration',false); 
            parse(p,varargin{:});       
            P.n = p.Results.refractive_index;
            P.method = p.Results.method;
            paraxial = p.Results.paraxial;
            P.gpu_acceleration = p.Results.gpu_acceleration;
            P.wavelength = wavelength_vacuum / P.n;
            P.k_prop = 2 * pi / P.wavelength;
            P.dist = dist;
            

            
            % Angular Spectrum (AS) and Band-Limited Angular Spectrum (BL-AS)
            if ismember(P.method, ["AS", "BL-AS"])
                if P.gpu_acceleration
                    grid_d2_fft_x2 = gpuArray((Grid.d2_fft_x .^ 2));
                    grid_d2_fft_y2 = gpuArray((Grid.d2_fft_y .^ 2));
                    P.mat = gpuArray(complex(zeros(size(grid_d2_fft_x2))));
                else
                    grid_d2_fft_x2 = Grid.d2_fft_x .^ 2;
                    grid_d2_fft_y2 = Grid.d2_fft_y .^ 2;
                    P.mat = zeros(size(grid_d2_fft_x2));
                    
                end
            elseif P.method == "FFT-DI"
                if P.gpu_acceleration
                    axis = gpuArray(complex(Grid.axis));
                else
                    axis = Grid.axis;
                end
            end
            
            if P.method == "AS"
               
                % 使用角谱理论
                condition = 1 / P.wavelength^2 - grid_d2_fft_x2 - grid_d2_fft_y2;
                if paraxial
                    P.mat = exp(1i * (P.k_prop * P.dist - pi * P.wavelength * (grid_d2_fft_x2 + grid_d2_fft_y2) * P.dist));
                else
                    P.mat(condition > 0) = exp(1i * 2 * pi * sqrt(max(condition, 0)) * P.dist);
                    P.mat(condition <= 0) = exp(-2 * pi * sqrt(abs(condition)) * P.dist);
                end
            elseif P.method == "BL-AS"
                % 使用带限角谱理论
                condition =  - grid_d2_fft_x2 - grid_d2_fft_y2 + 1 / P.wavelength^2;
                P.mat(condition > 0) = exp(1i * 2 * pi * sqrt(condition(condition > 0)) * P.dist);
                P.mat(condition <= 0) = 0;
                fft_limit_sq = 1 / ((2 * Grid.step_fft * P.dist)^2 + 1) / P.wavelength^2;
                condition1 = grid_d2_fft_x2 / fft_limit_sq + grid_d2_fft_y2 * P.wavelength^2 - 1;
                condition2 = grid_d2_fft_x2 * P.wavelength^2 + grid_d2_fft_y2 / fft_limit_sq - 1;
                P.mat((condition1 > 0) | (condition2 > 0)) = 0;  
            elseif P.method == "FFT-DI"
                % FFT Direct Integration (FFT-DI) method
                dr_real = sqrt(Grid.step^2 + Grid.step^2);
                rmax = sqrt(max(axis(:)).^2 + max(axis(:)).^2);
                dr_ideal = sqrt(P.wavelength^2 + rmax^2 + 2 * P.wavelength * sqrt(rmax^2 + P.dist^2)) - rmax;
                P.quality = dr_ideal / dr_real;
                if P.quality >= 1
                    fprintf('Good result: factor %.2f\n', P.quality);
                else
                    fprintf('Needs denser sampling: factor %.2f\n', P.quality);
                    fprintf('Distance: %.2f\n', P.dist);
                end
                X_tmp_1D = linspace(-2 * axis(1), 2 * axis(1), 2 * Grid.num_points - 1);
                [X_tmp_2D, Y_tmp_2D] = meshgrid(X_tmp_1D, X_tmp_1D);
                tmp_r = sqrt(X_tmp_2D.^2 + Y_tmp_2D.^2 + P.dist^2);
                P.mat_DI = (1 / (2 * pi)) * exp(1i * P.k_prop * tmp_r) .* (P.dist ./ tmp_r.^2) .* (1 ./ tmp_r - 1i * P.k_prop);
                P.mat_DI = P.mat_DI * Grid.step^2;
            end

        end
    end
end
