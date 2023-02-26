classdef reverberant_renderer < handle
    %REVERBERANT_RENDERER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fs
        c
        reverberant_source
        secondary_source_distribution
        output_signal
        delay_mx
        amp_mx
        delay_line
        prefilter
        antialiasing
      %  antialiasing_filters
    end
    
    methods
        function obj = reverberant_renderer(reverberant_source, SSD, fs, antialiasing)
            
            obj.antialiasing = antialiasing;
            obj.fs = fs;
            obj.c = 343.1;
            obj.reverberant_source = reverberant_source;
            obj.secondary_source_distribution = SSD;
            
            h_pre = obj.get_dim_mismatch_prefilter;
            obj.prefilter = OLS_convolver(h_pre, length(obj.reverberant_source.source_signal.time_series));
            for n = 1 : length(SSD)
                obj.output_signal{n} = signal;
            end
            obj = update_renderer( obj );

            N_blocks = length(obj.reverberant_source.source_signal.time_series);
            N = ceil(max(max(obj.delay_mx))/N_blocks);
            obj.delay_line = delay_line((N+1)*N_blocks,N_blocks, N);

        end
        
        function h = get_dim_mismatch_prefilter(obj)
            N = size(obj.reverberant_source.source_signal.time_series,1);
            w = [(0 : N/2 - 1)';(-N/2:-1)' ]/N*2*pi*obj.fs;
            H = sqrt(1i*w/obj.c);
            H(end/2+1) = real(H(end/2+1));
            h = fftshift(ifft(H));
            h = h(round(end/2)-50+1:round(end/2)+50).*hann(100);
        end
        
        function obj = update_renderer( obj )
            x0 =          cell2mat(cellfun( @(x) x.position, obj.secondary_source_distribution, 'UniformOutput', false)');
            n0 =          cell2mat(cellfun(  @(x) x.orientation, obj.secondary_source_distribution, 'UniformOutput', false)');
            center = mean(cell2mat(cellfun( @(x) x.position, obj.secondary_source_distribution, 'UniformOutput', false)'),1);
            rG = mean(sqrt(sum(bsxfun(@minus,cell2mat(cellfun( @(x) x.position, obj.secondary_source_distribution, 'UniformOutput', false)'),center).^2,2)));

            dl = mean(sqrt(sum((x0-circshift(x0,1)).^2,2)));
            
            xs = [ obj.reverberant_source.position cell2mat(cellfun( @(x) x.position, obj.reverberant_source.image_sources , 'UniformOutput', false))]';
            R_mx = zeros(size(xs,1),size(x0,1));
            
            for n = 1 : length(obj.secondary_source_distribution)
                k = bsxfun( @plus, x0(n,:), -xs);
                R_mx(:,n) = sqrt(sum(k.^2,2));
                k_P(:,:,n) = bsxfun( @times, k, 1./R_mx(:,n) );
                k_n(:,n) = squeeze(k_P(:,:,n))*n0(n,:)';
            end
            foc_ixs = find(all(k_n<=0,2));
            % Create delay matrix
            delay_mx = R_mx/obj.c*obj.fs;
            delay_mx(foc_ixs,:) = -delay_mx(foc_ixs,:);
            obj.delay_mx = ceil(delay_mx + ceil(rG/obj.c*obj.fs));

            % Create amplitude matrix
            obj.amp_mx = ((k_n).*double(k_n>=0))/sqrt(2*pi).*sqrt(R_mx.*rG./(R_mx + rG))./R_mx.*dl;
            % Treat focused sources as well
            xs = [ obj.reverberant_source.position cell2mat(cellfun( @(x) x.position, obj.reverberant_source.image_sources , 'UniformOutput', false))]';
            xs = xs(foc_ixs,:);
            k_P0 = bsxfun(@minus, xs,center);
            k_P0 = bsxfun(@times, k_P0, 1./sqrt(sum( k_P0.^2,2 )) );
            k_Pp = k_P(foc_ixs,:,:);
            w_0 = reshape(sum(k_Pp.*repmat(k_P0,1,1,size(k_Pp,3)),2),size(k_Pp,1),size(k_Pp,3));
            win = 0.5*(1-cos(2*pi*bsxfun(@times, w_0, 1./sum(w_0>0,2)) ));
            win = (bsxfun(@times, win, (1./max(win,[],2)) ).*(w_0>0)).^(1/4);
            
            obj.amp_mx(foc_ixs,:) = -win.*k_n(foc_ixs,:)/sqrt(2*pi).*...
                sqrt(R_mx(foc_ixs,:).*rG./(R_mx(foc_ixs,:) + rG))...
                ./R_mx(foc_ixs,:).*dl;
            obj.delay_mx(obj.amp_mx<=0) = nan;
        end
        
        
        function render(obj)
            obj.delay_line.write( obj.prefilter.convolve( obj.reverberant_source.source_signal.time_series) );
            for n = 1 : length(obj.secondary_source_distribution)
                obj.output_signal{n}.reset_signal;
                ixs = find(~isnan(obj.delay_mx(:,n)));
                for m = 1 : length(ixs)
                    obj.output_signal{n}.add_signals( obj.amp_mx(ixs(m),n)*...
                        obj.delay_line.read( obj.delay_mx(ixs(m),n),size(obj.reverberant_source.source_signal.time_series,1) ) );
                end
            end
        end
        
        
    end
end

