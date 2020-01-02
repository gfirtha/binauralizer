classdef freq_domain_delay_line < handle
    % Implements a frequency domain delay line (FDL) for an OLA or OLS
    % based partitioned convolver
    properties ( SetAccess = private )
        td_buffer
        td_ixs
        fdl
        fd_ixs
        N_fft
        N_part
        tap
    end
    properties ( SetAccess = private )
        write_ptr
    end
    methods
        function obj = freq_domain_delay_line( delay_size, partition_size, N_fft, N_blocks )
            obj.N_fft                = N_fft;
            [obj.td_ixs, obj.N_part] = obj.get_uniform_partitions( delay_size, partition_size );
            obj.td_buffer            = zeros(obj.td_ixs(end),1);
            obj.fd_ixs               = obj.get_uniform_partitions( N_fft/partition_size*delay_size, N_fft );
            obj.fdl                  = zeros(obj.fd_ixs(end),1);
            obj.write_ptr            = 1;
            obj.tap                  = cell(1,N_blocks);
        end
        
        function [block_ixs, N_uniform] = get_uniform_partitions(obj, delay_size, partition_size)
            N_uniform = ceil( delay_size / partition_size );
            partition_sizes = repmat( partition_size, N_uniform,1 );
            block_ixs = [ cumsum(partition_sizes)-partition_size+1, cumsum(partition_sizes)];
        end
        
        function output = read_td(obj, delay)
            output = obj.td_buffer( obj.td_ixs( ( 1 + mod( obj.write_ptr  - 1 -delay - 1, obj.N_part ) ) ,1 ):...
                obj.td_ixs( ( 1 + mod( obj.write_ptr  - 1 -delay - 1, obj.N_part ) ) ,2 ) );
        end
        function output = read_fd(obj, delay)
            output = obj.fdl( obj.fd_ixs( ( 1 + mod( obj.write_ptr  -2 -delay, obj.N_part ) ) ,1 ):...
                obj.fd_ixs( ( 1 + mod( obj.write_ptr  - 2 -delay, obj.N_part ) ) ,2 ) );
        end
        
        function obj = write(obj,input)
            % Some explanation:
            % 1 + mod( x , obj.N_part ): circular pointer, with 0 pointing
            % at the first elem
            
            obj.td_buffer(obj.td_ixs(obj.write_ptr,1):obj.td_ixs(obj.write_ptr,2)) = input;
            
            idx_prev_st = obj.td_ixs( ( 1 + mod( obj.write_ptr -2, obj.N_part ) ) ,1 );
            idx_prev_st = max( idx_prev_st, idx_prev_st + 2*length(input) - obj.N_fft );
            idx_prev_en = obj.td_ixs( ( 1 + mod( obj.write_ptr -2, obj.N_part ) ) ,2 );
            idx_st = obj.td_ixs( ( 1 + mod( obj.write_ptr -1, obj.N_part ) ) ,1 );
            idx_en = obj.td_ixs( ( 1 + mod( obj.write_ptr -1, obj.N_part ) ) ,2 );
            
            obj.fdl(obj.fd_ixs(obj.write_ptr,1):obj.fd_ixs(obj.write_ptr,2)) =  ...
               fft([obj.td_buffer( idx_prev_st : idx_prev_en); ...
                    obj.td_buffer( idx_st      :      idx_en)], obj.N_fft);
%             obj.fdl(obj.fd_ixs(obj.write_ptr,1):obj.fd_ixs(obj.write_ptr,2)) =  ...
%                 fft(obj.td_buffer(1 + mod( idx_en - 1 + (-min(obj.N_fft,2*length(input))+1:0), ...
%                 length(obj.td_buffer) )), obj.N_fft);
            for n = 1 : obj.N_part
                obj.tap{n} = obj.fd_ixs( 1 + mod( obj.write_ptr -n, obj.N_part ), : );
            end
            obj.write_ptr = 1 + mod( obj.write_ptr , obj.N_part );
        end
    end
end