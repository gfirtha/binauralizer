classdef delay_line < handle
    % Implements a frequency domain delay line (FDL) for an OLA or OLS
    % based partitioned convolver
    properties ( SetAccess = private )
        td_buffer
        td_ixs
        N_part
    end
    properties ( SetAccess = private )
        write_ptr
    end
    methods
        function obj = delay_line( delay_size, partition_size, ~ )
            [obj.td_ixs, obj.N_part] = obj.get_uniform_partitions( delay_size, partition_size );
            obj.td_buffer            = zeros(obj.td_ixs(end),1);
            obj.write_ptr            = 1;
        end
        
        function [block_ixs, N_uniform] = get_uniform_partitions(~, delay_size, partition_size)
            N_uniform = ceil( delay_size / partition_size );
            partition_sizes = repmat( partition_size, N_uniform,1 );
            block_ixs = [ cumsum(partition_sizes)-partition_size+1, cumsum(partition_sizes)];
        end
        
        function output = read(obj, delay, block_size)
            start_idx = 1 + mod(  obj.td_ixs( 1 + mod( obj.write_ptr - 2,obj.N_part ) ) - round(delay) -1 , length(obj.td_buffer) );
            end_idx = 1 + mod( start_idx - 1 + block_size -1 , length(obj.td_buffer) );
            if start_idx > end_idx
               output = [obj.td_buffer( start_idx:end ); obj.td_buffer( 1:end_idx)];
            else
               output = obj.td_buffer( start_idx : end_idx );
            end
        end
        
        function obj = write(obj,input)
            % Some explanation:
            % 1 + mod( x , obj.N_part ): circular pointer, with 0 pointing
            % at the first elem
            
            obj.td_buffer(obj.td_ixs(obj.write_ptr,1):obj.td_ixs(obj.write_ptr,2)) = input;
            obj.write_ptr = 1 + mod( obj.write_ptr , obj.N_part );
        end
    end
end