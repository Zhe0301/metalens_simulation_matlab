function  plane_wave(S, alpha, beta)
    % 平面波函数
    % alpha: 是平面波x方向余弦，弧度制
    % beta: 是平面波y方向余弦，弧度制
    S.k_prop = 2 * pi / S.wavelength_vacuum;
    if alpha == pi / 2 && beta == pi / 2
        S.phase = zeros(size(S.Grid.d2_x));
    else
        S.phase = S.k_prop * (S.Grid.d2_x * cos(alpha) + S.Grid.d2_y * cos(beta));
    end
    S.complex_amplitude = S.amplitude .* exp(1i * S.phase);
end