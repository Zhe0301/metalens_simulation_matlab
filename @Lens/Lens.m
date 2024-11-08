classdef Lens < handle
    properties
        Grid
        complex_amplitude_t
        phase
        mask
        amplitude
    end
    
    methods
        function L = Lens(Grid)
            L.Grid = Grid;
            L.complex_amplitude_t = [];
            L.phase = [];
            L.mask = [];
            L.amplitude = [];
        end

        moire_quadratic(L, r_max, a, phi,varargin)
        binary2(L, r_max, a, varargin)
        binary2_d(L, r_max, a, varargin)
        ideal_lens(L, r_max, focal_length, wavelength_vacuum, varargin)
        hole(L, r_max, varargin)
        square_hole(L, length_x, length_y, varargin)
        plot_phase(L, varargin)
        plot_intensity(L, varargin)
    end
end
