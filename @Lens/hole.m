function hole(L, r_max, varargin)
    % 圆孔
    % r_max: 孔的半径，单位mm;
    p = inputParser;    
    addParameter(p,'t',1); 
    parse(p,varargin{:});       
    t = p.Results.t;
    
    mask_index = L.Grid.d2_r <= r_max;
    L.mask = zeros(size(L.Grid.d2_r));
    L.mask(mask_index) = 1;
    L.amplitude = ones(size(L.Grid.d2_r)) .* L.mask;
    L.phase = zeros(size(L.Grid.d2_r)) .* L.mask;
    L.complex_amplitude_t = L.amplitude .* exp(1i * L.phase) * t .* L.mask;
end
