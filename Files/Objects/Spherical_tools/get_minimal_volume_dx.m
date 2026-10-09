        function [r,kr0] = get_minimal_volume_dx(obj,k0)
            k0 = k0./sqrt(sum(k0.^2,1));
            kr = k0'-sum(obj.voronoi_cells.n_in.*k0',2).*obj.voronoi_cells.n_in;
            kr0 = sqrt(sum(kr.^2,2));
            kr = kr./kr0;
            for n = 1 : size(kr,1)
                wtr(:,n) = squeeze(obj.min_vol_ellipses.Tmx(:,:,n))*(kr(n,:))';
            end
            r = ( prod(obj.min_vol_ellipses.R,1)./sqrt(sum((obj.min_vol_ellipses.R.*flipud(wtr)).^2,1)) )';
            r(isnan(r)) = 1e-10;
        end
