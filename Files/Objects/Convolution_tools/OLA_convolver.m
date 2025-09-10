classdef OLA_convolver < base_fft_convolver
    %OLA_CONVOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = OLA_convolver(coefficients)
            obj.N_filt = size( coefficients, 1 );
            obj.N_fft  = size( coefficients, 1 );
            if size( coefficients, 2 ) == 1
                obj.coefficients = fft( coefficients, obj.N_fft, 1);
            else
                obj.coefficients = fft( coefficients(:,1) + 1i*coefficients(:,2) , obj.N_fft, 1);
            end
            obj.state_buffer = zeros( obj.N_filt - 1 , 1 );
        end
        
        function obj = update_coefficients(obj,coefficients)
            obj.N_filt = size( coefficients, 1 );
            if size( coefficients, 2 ) == 1
                obj.coefficients = fft( coefficients, obj.N_fft, 1);
            else
                obj.coefficients = fft( coefficients(:,1) + 1i*coefficients(:,2) , obj.N_fft, 1);
            end
        end
        
        function output = convolve(obj, input)
            obj.check_fft_size( size( input, 1) );
            input_fft =  fft( input, obj.N_fft );
            out_full  = ifft( input_fft.*obj.coefficients );
            output    = out_full(1:length(input)) + ...
                [obj.state_buffer; zeros( length(input) - length(obj.state_buffer), 1)];
            obj.state_buffer = out_full( length(input) + 1 : length(input) + obj.N_filt - 1 );
            if ~isreal(output)
                output = [real(output), imag(output)];
            end
        end
        
        function out_full = partitioned_convolve(obj,input_fft,prev_state_buffer)
            out_full  = ifft( input_fft.*obj.coefficients );
        end
    end
end

