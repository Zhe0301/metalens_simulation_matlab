function [mtf_x, mtf_y, psf] = calculate_MTF(psf, step, varargin)
% 计算 MTF
% 参数:
% psf: 点扩散函数
% step: 步长 (mm)
% 可选参数:
% zero_coef: 置零系数，小于 max(psf)/zero_coef 的值视为 0, 默认为1e5
% gpu_acceleration: 是否使用 GPU 加速计算, 默认为false
% 返回值:
% mtf_x: x 方向的调制传递函数
% mtf_y: y 方向的调制传递函数
% psf: 置零处理后的 psf
    p = inputParser;            
    addParameter(p,'zero_coef',1e5);      
    addParameter(p,'gpu_acceleration',false); 
    parse(p,varargin{:}); 
    zero_coef = p.Results.zero_coef;
    gpu_acceleration = p.Results.gpu_acceleration;
    if gpu_acceleration
        psf = gpuArray(psf); % 转换为 GPU 数组
    end
    % 将小于 max(psf) / zero_coef 的值置为零
    max_psf = max(psf(:));
    psf(psf < max_psf / zero_coef) = 0;

    % 计算 x 和 y 方向的线扩散函数 (LSF)
    lsf_x = trapz(psf, 1) * step; % 对 y 积分
    lsf_y = trapz(psf, 2) * step; % 对 x 积分

    % 计算 MTF（调制传递函数）
    mtf_x = abs(fftshift(fft(lsf_x)));
    mtf_y = abs(fftshift(fft(lsf_y)));

    % 只保留正频部分并归一化
    N = numel(mtf_x);
    mtf_x = mtf_x(N/2+1:end);
    mtf_x = mtf_x / max(mtf_x);
    
    mtf_y = mtf_y(N/2+1:end);
    mtf_y = mtf_y / max(mtf_y);

    % 如果使用了 GPU 加速，则将结果转换为普通数组
    if gpu_acceleration
        mtf_x = gather(mtf_x);
        mtf_y = gather(mtf_y);
        psf = gather(psf);
    end
end
