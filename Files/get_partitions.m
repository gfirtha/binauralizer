function ix = get_partitions( N, block_size, mode )
        N_part = ceil( N/block_size );
        n = (0:N_part-1)';
        ix = [ n*block_size + 1, (n+1)*block_size ];
        ix(ix>N) = N;
            
end

