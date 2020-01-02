classdef base_fft_convolver < base_convolver

    properties (SetAccess = protected)
        N_fft
    end
    
    methods 
          function obj = check_fft_size(obj, input_size)
            if obj.N_fft < ( obj.N_filt + input_size - 1 )
               obj.N_fft = 2^nextpow2( obj.N_filt + input_size - 1 );
               obj.coefficients = fft( ifft(obj.coefficients), obj.N_fft, 1 );
            end
          end
    end
        
end

