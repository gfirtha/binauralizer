classdef OLS_convolver < base_fft_convolver
    
    properties
        delay
        block_size
        part_ix
        mode
    end
    
    methods
        function obj = OLS_convolver(coefficients, block_size)
            obj.N_filt = min( size( coefficients, 1 ), block_size );
            obj.block_size = block_size;
            obj.N_fft =  2^nextpow2( min(block_size + size( coefficients, 1 ), 2*block_size) - 1 );
            if size(coefficients,1) <= block_size
                obj.mode = 'short_convolution';
                obj.state_buffer = zeros( obj.N_filt - 1 , 1 );
            else
                obj.mode = 'uniform_partitions';
                N_blocks = ceil(size(coefficients,1)/block_size);
                if N_blocks == 1
                    obj.delay = freq_domain_delay_line( (N_blocks + 1)*block_size, ...
                        block_size, obj.N_fft, N_blocks );
                else
                    obj.delay = freq_domain_delay_line( N_blocks*block_size, ...
                        block_size, obj.N_fft, N_blocks );
                end
            end
            obj.part_ix = get_partitions(size(coefficients,1), block_size, 'uniform');
            for n = 1 : size(obj.part_ix,1)
                if size( coefficients, 2 ) == 1
                    obj.coefficients{n} = fft( coefficients(obj.part_ix(n,1):obj.part_ix(n,2)), obj.N_fft, 1);
                else
                    obj.coefficients{n} = fft( coefficients(obj.part_ix(n,1):obj.part_ix(n,2),1) +...
                        1i*coefficients(obj.part_ix(n,1):obj.part_ix(n,2),2) , obj.N_fft, 1);
                end
            end
        end
        
        function obj = update_coefficients(obj,coefficients)
            for n = 1 : size(obj.part_ix,1)
                if size( coefficients, 2 ) == 1
                    obj.coefficients{n} = fft( coefficients(obj.part_ix(n,1):obj.part_ix(n,2)), obj.N_fft, 1);
                else
                    obj.coefficients{n} = fft( coefficients(obj.part_ix(n,1):obj.part_ix(n,2),1) +...
                        1i*coefficients(obj.part_ix(n,1):obj.part_ix(n,2),2) , obj.N_fft, 1);
                end
            end
        end
        
        function output = convolve(obj, input)
            switch obj.mode
                case 'uniform_partitions'
                    obj.delay.write(input);
                    output_fft = 0;
                    for n = 1 : length(obj.coefficients)
                        output_fft = output_fft + obj.coefficients{n} .* obj.delay.fdl( obj.delay.tap{n}(1):obj.delay.tap{n}(2)  );
                    end
                    output_full = ifft( output_fft );
                    if obj.N_fft <= 2*length(input)
                        output = output_full( obj.N_fft - length(input) + 1  : end);
                    else
                        output = output_full( length(input) + 1  : 2*length(input));
                    end
                case 'short_convolution'
                    input_fft =  fft( [ obj.state_buffer; input], obj.N_fft );
                    out_full  = ifft( input_fft.*obj.coefficients{1} );
                    obj.state_buffer = input(end - obj.N_filt + 2 : end);
                    output    = out_full( obj.N_filt : obj.N_filt + size(input,1) -1 );
            end
            if ~isreal(output)
                output = [real(output), imag(output)];
            end
        end
        
    end
end

