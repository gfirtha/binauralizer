classdef base_convolver < handle

    properties ( SetAccess = protected )
        coefficients
        state_buffer
        N_filt
    end
    
    methods
        function obj = clear_coeffs(obj)
            obj.coefficients = zeros( size( obj.coefficients ) );
        end
        
        function obj = clear_buffer(obj)
            obj.state_buffer = zeros( obj.N_filt - 1 , 1 );
        end
    end
    
    methods (Abstract)
        update_coefficients(obj, coefficients)
        output = convolve(obj, input)
    end
end

