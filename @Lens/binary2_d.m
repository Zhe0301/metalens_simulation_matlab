function binary2_d(L, r_max, a, varargin)
    % 相位离散型的二元面2模拟: 使用单元结构的真实相位进行仿真
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 参考zemax二元面
    % r_max：最大半径,单位 mm;
    % a: 多项式系数数组;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 可选参数
    % m: 衍射阶数，默认值为1;
    % r_0: 归一化半径,单位 mm, 默认值为1;
    % unit_phase(N维数组): 离散单元的相位值
    % unit_t(N维数组):离散单元的能量透过率
    % boundaries(2*N矩阵):每个离散相位值代替的相位范围
    p = inputParser;            
    addParameter(p,'m',1);      
    addParameter(p,'r_0',1); 
    addParameter(p,'unit_phase',(0:7) * pi/4 + pi/8); 
    addParameter(p,'unit_t',ones(1, 8)); 
    addParameter(p,'boundaries',[0, 1/8; 1/8, 2/8; 2/8, 3/8; 3/8, 4/8; ...
        4/8, 5/8; 5/8, 6/8; 6/8, 7/8; 7/8, 1] * 2 * pi); 
    parse(p,varargin{:}); 
    tic;
    m = p.Results.m;
    r_0 = p.Results.r_0;
    unit_phase = p.Results.unit_phase;
    unit_t = p.Results.unit_t;
    boundaries = p.Results.boundaries;
    
    mask_index = L.Grid.d2_r <= r_max;
    L.mask = zeros(size(L.Grid.d2_r));
    L.mask(mask_index) = 1;
    L.phase = zeros(size(L.Grid.d2_r));
    L.amplitude = ones(size(L.Grid.d2_r)) .* L.mask;

    for i = 1:length(a)
        L.phase = L.phase + m * a(i) * (L.Grid.d2_r / r_0) .^ (2 * i);
    end
    
    L.phase = mod(L.phase, 2 * pi);
    
    if isempty(unit_phase) || isempty(unit_t) || isempty(boundaries)
        disp("Error: Need proper discrete values for phase and transmittance.");
        L.phase = L.phase .* L.mask;
    else
        for i = 1:size(boundaries, 1)
            mask = (L.phase >= boundaries(i, 1)) & (L.phase < boundaries(i, 2));
            L.phase(mask) = unit_phase(i);
            L.amplitude(mask) = sqrt(unit_t(i));
        end
        L.phase = L.phase .* L.mask;
    end

    L.complex_amplitude_t = L.amplitude .* exp(1i * L.phase) .* L.mask;
    elapsedTime = toc;
    fprintf("binary2_d initialization complete: %.2fs\n",elapsedTime);
end
