classdef OLS_convolver < base_fft_convolver
    
    properties
        delay
        block_size
        part_ix
        mode
        prefilter
        unfiltered_coefficients
        output_mode
    end
    
    methods
        function obj = OLS_convolver(varargin)
            coefficients = varargin{1};
            obj.block_size = varargin{2};
            if nargin > 2
                obj.output_mode = varargin{3};
            else
                obj.output_mode = 'time_signal';
            end
            
            obj.N_filt = min( size( varargin{1}, 1 ), varargin{2} );
            obj.N_fft =  2^nextpow2( min(obj.block_size + size( coefficients, 1 ), 2*obj.block_size) - 1 );
            if size(coefficients,1) <= obj.block_size
                obj.mode = 'short_convolution';
                obj.state_buffer = zeros( obj.N_filt - 1 , 1 );
            else
                obj.mode = 'uniform_partitions';
                N_blocks = ceil(size(coefficients,1)/obj.block_size);
                if N_blocks == 1
                    obj.delay = freq_domain_delay_line( (N_blocks + 1)*obj.block_size, ...
                        obj.block_size, obj.N_fft, N_blocks );
                else
                    obj.delay = freq_domain_delay_line( N_blocks*obj.block_size, ...
                        obj.block_size, obj.N_fft, N_blocks );
                end
            end
            obj.part_ix = get_partitions(size(coefficients,1), obj.block_size, 'uniform');
            for n = 1 : size(obj.part_ix,1)
                if size( coefficients, 2 ) == 1
                    obj.coefficients{n} = fft( coefficients(obj.part_ix(n,1):obj.part_ix(n,2)), obj.N_fft, 1);
                else
                    obj.coefficients{n} = fft( coefficients(obj.part_ix(n,1):obj.part_ix(n,2),1) +...
                        1i*coefficients(obj.part_ix(n,1):obj.part_ix(n,2),2) , obj.N_fft, 1);
                end
            end
            obj.unfiltered_coefficients = obj.coefficients;
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
            obj.unfiltered_coefficients = obj.coefficients;
            if ~isempty(obj.prefilter) 
                obj.prefilter_coefficients;
            end
        end
        
        function obj = update_prefilter(obj,prefilter,input_mode)
            switch input_mode
                case 'frequency_domain'
                    obj.prefilter = prefilter;
                case 'time_domain'
                    obj.prefilter = fft(prefilter, obj.N_fft );
            end
            obj.prefilter_coefficients;
        end
        function obj = prefilter_coefficients(obj)
            % TODO: extend prefilter for UPOLS convolution
            if ~all(obj.prefilter==1)
                for n = 1 : size(obj.part_ix,1)
                    if size( obj.coefficients, 2 ) == 1
                        obj.coefficients{n} = obj.prefilter.*obj.unfiltered_coefficients{1};
                    else
                        obj.coefficients{n} = bsxfun( @times, obj.prefilter, obj.unfiltered_coefficients(obj.part_ix(n,:):obj.part_ix(n,:)) ) ;
                    end
                end
%            else 
            end
        end
        
        function varargout = convolve(obj, input)
            switch obj.mode
                case 'uniform_partitions'
                    obj.delay.write(input);
                    output_fft = 0;
                    for n = 1 : length(obj.coefficients)
                        output_fft = output_fft + obj.coefficients{n} .* obj.delay.fdl( obj.delay.tap{n}(1):obj.delay.tap{n}(2)  );
                    end
                    switch obj.output_mode
                        case 'time_signal'
                            output_full = ifft( output_fft );
                            if obj.N_fft <= 2*length(input)
                                output = output_full( obj.N_fft - length(input) + 1  : end);
                            else
                                output = output_full( length(input) + 1  : 2*length(input));
                            end
                            if ~isreal(output)
                                output = [real(output), imag(output)];
                            end
                            varargout{1} = output;
                        case 'spectrum'
                            
                            if obj.N_fft <= 2*length(input)
                                ix = [obj.N_fft - length(input) + 1 , length(output_fft)];
                            else
                                ix = [length(input) + 1  , 2*length(input)];
                            end
                            varargout{1} = output_fft;
                            varargout{2} = ix;
                            
                    end
                case 'short_convolution'
                    input_fft =  fft( [ obj.state_buffer; input], obj.N_fft );
                    obj.state_buffer = input(end - obj.N_filt + 2 : end);
                    switch obj.output_mode
                        case 'time_signal'
                            out_full  = ifft( input_fft.*obj.coefficients{1} );
                            output    = out_full( obj.N_filt : obj.N_filt + size(input,1) -1 );
                            if ~isreal(output)
                                output = [real(output), imag(output)];
                            end
                            varargout{1} = output;
                        case 'spectrum'
                            output = input_fft.*obj.coefficients{1};
                            ix = [obj.N_filt, (obj.N_filt + size(input,1) -1)];
                            varargout{1} = output;
                            varargout{2} = ix;
                    end
            end
        end
    end
end
