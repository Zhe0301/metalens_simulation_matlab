classdef Source < handle
    properties
        Grid
        phase
        k_prop
        complex_amplitude
        wavelength_vacuum
        amplitude
    end
    
    methods
        function S = Source(Grid, wavelength_vacuum, amplitude)
            S.Grid = Grid;
            S.phase = [];
            S.k_prop = [];
            S.complex_amplitude = [];
            S.wavelength_vacuum = wavelength_vacuum;
            S.amplitude = amplitude * ones(size(Grid.d2_x));
        end
        plane_wave(S, alpha, beta)
        plot_intensity(S, varargin)
        plot_phase(S, varargin)
        
    end
end
