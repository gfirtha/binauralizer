        function [quadrature_points,output,amp,phase] = interpolate2quadrature(input,delauney_trias,x0)
            if ~isreal(input)
                A_in = abs(input);
                phase_in = unwrap(angle(input),[],3);
            end
            quadrature_points = get_adaptive_quadrature_points(delauney_trias,x0);
            sizein = size(input);
            output = zeros([size(quadrature_points.quad_pos,1),sizein(2:end)]);
            amp = zeros([size(quadrature_points.quad_pos,1),sizein(2:end)]);
            phase = zeros([size(quadrature_points.quad_pos,1),sizein(2:end)]);
            for n = 1 : size(quadrature_points.quad_pos,1)
                if  isreal(input)
                    output(n,:,:) =  sum(quadrature_points.interp_w(n,:)'.*input(quadrature_points.interp_pos(n,:)',:,:),1);
                else
                    A_interp = sum(quadrature_points.interp_w(n,:)'.*A_in(quadrature_points.interp_pos(n,:)',:,:),1);
                    phase_interp = sum(quadrature_points.interp_w(n,:)'.*phase_in(quadrature_points.interp_pos(n,:)',:,:),1);
                    output(n,:,:) =  A_interp.*exp(1i*phase_interp);
                    amp(n,:,:) =  A_interp;
                    phase(n,:,:) =  phase_interp;
                end
            end
        end
