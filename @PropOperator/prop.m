function e_out = prop(P, complex_amplitude)
    if P.gpu_acceleration
        complex_amplitude = gpuArray(complex(complex_amplitude));
    end
    
    if ismember(P.method, ["AS", "BL-AS"])
        wave_fft = fftshift(fft2(complex_amplitude));
        wave_prop = wave_fft .* P.mat;
        P.e_out = ifft2(ifftshift(wave_prop));
    elseif P.method == "FFT-DI"
        N = size(complex_amplitude, 1);
        Mat_U = zeros(2 * N - 1, 2 * N - 1, 'like', complex_amplitude);
        Mat_U(1:N, 1:N) = complex_amplitude;
        S = ifft2(fft2(Mat_U) .* fft2(P.mat_DI));
        P.e_out = S(N:end, N:end);
    end
    
    if P.gpu_acceleration
        P.e_out = gather(P.e_out);
    end
    e_out = P.e_out;
end