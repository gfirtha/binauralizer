        function plot_ellipse(obj,n)
            u = linspace(0, 2*pi, 90);
            X_ell = pinv(obj.min_vol_ellipses.Tmx(:,:,n)) * [ obj.min_vol_ellipses.R(1,n)*cos(u);...
                obj.min_vol_ellipses.R(2,n)*sin(u) ] + obj.min_vol_ellipses.center(:,n);
            hold on
            plot3(X_ell(1,:)',X_ell(2,:)',X_ell(3,:)','k','LineWidth',1.5);
        end
