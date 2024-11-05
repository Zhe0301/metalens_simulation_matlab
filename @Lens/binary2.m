function binary2(L, r_max, a, varargin)
    % 二元面2模拟: 使用径向坐标的偶次多项式描述相位面
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 参考zemax二元面
    % r_max：最大半径,单位 mm;
    % a: 多项式系数数组;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 可选参数
    % m: 衍射阶数，默认值为1;
    % r_0: 归一化半径,单位 mm, 默认值为1;
    % t: 能量透过率，默认为1
    p = inputParser;            
    addParameter(p,'m',1);      
    addParameter(p,'r_0',1); 
    addParameter(p,'t',1); 
    parse(p,varargin{:});       
    m = p.Results.m;
    r_0 = p.Results.r_0;
    t = p.Results.t;
    
    mask_index = L.Grid.d2_r <= r_max;
    L.mask = zeros(size(L.Grid.d2_r));
    L.mask(mask_index) = 1;
    L.phase = zeros(size(L.Grid.d2_r));
    
    for i = 1:length(a)
        L.phase = L.phase + m * a(i) * (L.Grid.d2_r / r_0) .^ (2 * i);
    end
    
    L.phase = L.phase .* L.mask;
    L.amplitude = ones(size(L.Grid.d2_r)) .* L.mask;
    L.complex_amplitude_t = L.amplitude .* exp(1i * L.phase) * t .* L.mask;
end
