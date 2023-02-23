function [hrtf_resp] = get_hrtfs( x_s, x_r, o_r, hrtf_db, hrtf_spectra )

interp_method = 'spline';
x_extrap = bsxfun( @minus, x_s, x_r );
R_extrap = sqrt(sum(x_extrap.^2,2));
R_measurement = mean(hrtf_db.SourcePosition(:,3));
v_sr = bsxfun(@times, x_extrap, 1./R_extrap);
[theta_1,~] = cart2pol(v_sr(:,1), v_sr(:,2));
[theta_2,~] = cart2pol(o_r(1),o_r(2));

hrtf_resp = zeros(2, size(x_s,1),size(hrtf_db.Data.IR,3));
for n = 1 : size(x_extrap,1)
    theta_interp = mod( theta_1(n) - theta_2, 2*pi);
    hrtf_interp = interpolate_hrtf( hrtf_spectra.spectrum, hrtf_spectra.theta, theta_interp, interp_method);
%    [theta_ext, hrtf_ext] = append_hrtf(hrtf_spectra.theta, hrtf_spectra.spectrum, hrtf_interp, theta_interp );

%    hrtf_extrap = extrapolate_hrtf( hrtf_ext, hrtf_db.Data.SamplingRate ,x_extrap(n,:),...
%        R_measurement*[cos(theta_ext) sin(theta_ext)]);

    hrtf_resp(:,n,:) = ifft( hrtf_interp, [], 2,'symmetric' );
end
hrtf_resp = squeeze( hrtf_resp );

end

