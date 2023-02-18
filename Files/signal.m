classdef signal < handle
    %SIGNAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        time_series
        spectrum
        opt_indexes
        fs
    end
    
    methods
        function obj = signal(varargin)
            if nargin > 0
                obj.time_series = varargin{1};
            end
            if nargin > 1
                obj.fs = varargin{2};
            end

        end
        
        function obj = set_signal(obj, varargin)
            obj.time_series = varargin{1};
            if length(varargin)>1
                obj.fs = varargin{2};
            end
        end
        function obj = set_spectrum(obj, varargin)
            if nargin == 1
                obj.spectrum = varargin{1};
            else
                obj.spectrum = varargin{1};
                obj.opt_indexes = varargin{2};
            end
        end
        function output = get_signal(obj)
            if ~isempty(obj.time_series)
                output = obj.time_series;
            else
                if isempty(obj.opt_indexes)
                    output = ifft(obj.spectrum);
                else
                    out_temp = ifft(obj.spectrum);
                    output = out_temp(obj.opt_indexes(1):obj.opt_indexes(2));
                end
            end
            if ~isreal(output)
                output = [real(output), imag(output)];
            end
        end
        function output = get_spectrum(obj)
            if ~isempty(obj.spectrum)
                output = obj.spectrum;
            else
                output = fft(obj.time_series);
            end
        end
        
        function obj = add_signals(obj, input)
            
            if isempty(obj.time_series)
                obj.time_series = input.time_series;
            else
                obj.time_series = obj.time_series + input.time_series;
            end
                
        end
        
        function obj = add_spectra(obj, input)
            if isempty(obj.opt_indexes)
                obj.opt_indexes = input.opt_indexes;
            end
            if isempty(obj.spectrum)
                obj.spectrum = input.spectrum;
            else
                obj.spectrum = obj.spectrum + input.spectrum;
            end
        end
    end
end

