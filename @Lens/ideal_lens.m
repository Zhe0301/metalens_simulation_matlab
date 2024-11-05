function ideal_lens(L, r_max, focal_length, wavelength_vacuum, varargin)
    % 理想透镜函数
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % r_max: 镜片最大半径;
    % focal_length: 焦距;
    % wavelength_vacuum: 真空中波长
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 可选参数
    % t: 能量透过率  
    p = inputParser;    
    addParameter(p,'t',1); 
    parse(p,varargin{:});       
    t = p.Results.t;
    
    mask_index = L.Grid.d2_r <= r_max;
    L.mask = zeros(size(L.Grid.d2_r));
    L.mask(mask_index) = 1;
    L.phase = zeros(size(L.Grid.d2_r));
    
    L.phase = -(focal_length - sqrt(focal_length^2 + L.Grid.d2_r.^2)) * 2 * pi / wavelength_vacuum;
    L.phase = L.phase .* L.mask;
    L.amplitude = ones(size(L.Grid.d2_r)) .* L.mask;
    L.complex_amplitude_t = L.amplitude .* exp(-1i * L.phase) * t .* L.mask;
end