function square_hole(L, length_x, length_y, varargin)
    % 方孔
    % lenth_x: 孔的x长度;
    % lenth_y: 孔的y长度;
    p = inputParser;
    addParameter(p,'t',1); 
    parse(p,varargin{:});       
    t = p.Results.t;
    
    mask_index = (abs(L.Grid.d2_x) <= length_x / 2) & (abs(L.Grid.d2_y) <= length_y / 2);
    L.mask = zeros(size(L.Grid.d2_r));
    L.mask(mask_index) = 1;
    L.amplitude = ones(size(L.Grid.d2_r)) .* L.mask;
    L.phase = zeros(size(L.Grid.d2_r)) .* L.mask;
    L.complex_amplitude_t = L.amplitude .* exp(1i * L.phase) * t .* L.mask;
end