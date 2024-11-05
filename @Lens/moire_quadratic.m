function moire_quadratic(L, r_max, a, phi,varargin)
    %二次相位摩尔透镜模拟
    % 参考文献：Adjustable refractive power from diffractive moire elements (2008)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % r_max: 最大半径,单位 mm;
    % a: 相位系数，影响变焦范围.a的正负可区分第一面和第二面
    % phi: 旋转角度，弧度制，换算值0-2pi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 可选参数
    % round_off: 是否使用相位四舍五入近似，默认为true
    % f_offset_wavelength: 额外的焦距补偿项，是补偿焦距和波长的乘积。补偿焦距是旋转角度为零时的焦距。为0时无补偿。默认为0
    % t: 能量透过率，默认为1
    p = inputParser;            
    addParameter(p,'f_offset_wavelength',0);      % 设置变量名和默认参数
    addParameter(p,'round_off',true); 
    addParameter(p,'t',1); 
    parse(p,varargin{:});       % 对输入变量进行解析，如果检测到前面的变量被赋值，则更新变量取值
    f_offset_wavelength = p.Results.f_offset_wavelength;
    round_off = p.Results.round_off;
    t = p.Results.t;
    
    mask_index = L.Grid.d2_r <= r_max;
    L.mask = zeros(size(L.Grid.d2_r));
    L.mask(mask_index) = 1;
    L.phase = zeros(size(L.Grid.d2_r));
    differential_angle = L.Grid.d2_theta - phi;
    
    if round_off
        obj.phase = ceil(a * L.Grid.d2_square) .* differential_angle;
    else
        differential_angle(differential_angle <= -pi) = differential_angle(differential_angle <= -pi) + 2 * pi;
        L.phase = a * L.Grid.d2_square .* differential_angle;
    end
    
    if f_offset_wavelength ~= 0
        L.phase = L.phase + L.Grid.d2_square * pi / (2 * f_offset_wavelength);
    end
    L.phase = L.phase .* L.mask;
    L.amplitude = ones(size(L.Grid.d2_r)) .* L.mask;
    L.complex_amplitude_t = L.amplitude .* exp(1i * L.phase) * t .* L.mask;
    end
