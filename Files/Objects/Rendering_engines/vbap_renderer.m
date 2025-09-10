classdef vbap_renderer < base_renderer
    %VBAP_RENDERER Summary of this class goes here
    %   Detailed explanation goes here

    properties
        G_vec
    end

    methods
        function obj = vbap_renderer(virtual_source,loudspeaker_array)
            obj = obj@base_renderer(virtual_source,loudspeaker_array);
            obj.virtual_source = virtual_source;
%            obj.ssd.loudspeakers = loudspeaker_array;
            obj.update_renderer;
        end
        function info = get_renderer_info(obj)
            info = "The renderer implements simple Vector Based Amplitude Panning between loudspeaker pairs.";
        end

        function obj = update_renderer(obj, ~)
            obj.G_vec = zeros(length(obj.ssd.loudspeakers),1);
            L = cell2mat(cellfun( @(x) x.position', obj.ssd.loudspeakers,'UniformOutput',false));
            R = cell2mat(cellfun( @(x) sqrt(sum(x.position.^2,2)), obj.ssd.loudspeakers,'UniformOutput',false));
            L = bsxfun(@times, L, 1./R);
            % fi = cart2pol(L(1,:),L(2,:));
            % %            [~,i] = sort(abs( cart2pol(obj.virtual_source.positionobj.virtual_source.position(1),obj.virtual_source.position(2)) - fi ));
            % [val,i] = sort(L'*obj.virtual_source.position','descend');
            % l = L(:,i(1:2));
            % ind = i(1:2);
            % g = l\obj.virtual_source.position'/norm(obj.virtual_source.position);
            % obj.G_vec(ind) = g/norm(g);
            ls_pairs = [(1:length(obj.ssd.loudspeakers));
                circshift((1:length(obj.ssd.loudspeakers)),1)];
            xvs = obj.virtual_source.position*norm(obj.virtual_source.position);
            for n = 1 : size(ls_pairs,2)
                g(:,n) = ( L(:,ls_pairs(:,n)) \ xvs' );
            end
            ixs = find(sum( (g>=0), 1)==2);
            dl =  sqrt(sum((obj.ssd.x0-circshift(obj.ssd.x0,1)).^2,2));

            [~,ix] = min(dl(ixs));
            ix = ixs(ix);
        
            R = sqrt( sum( (obj.ssd.x0).^2,2 ) );
            obj.G_vec(ls_pairs(:,ix)) = g(:,ix)/norm(g(:,ix));

        end

        function render(obj)
            for n = 1 : length(obj.ssd.loudspeakers)
                obj.output_signal{n}.set_signal( obj.G_vec(n)*obj.virtual_source.source_signal.time_series );
            end
            obj.add_output_to_ssd_signal;
        end

    end
end