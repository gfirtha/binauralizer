function [hrtf_interp] = interpolate_hrtf(hrtf_meas, theta_measurement,  theta_interp, mode)
%%
[~,ix] = min( abs( theta_measurement - theta_interp ) );
switch mode
    case 'nearest'
    % Nearest HRTF interpoaltion
    hrtf_interp = hrtf_meas(ix, :, :);

    case 'spline'
    % Direct cubic interpoaltion in frequency domain
    hrtf_l = squeeze(hrtf_meas(:,1,:));
    hrtf_r = squeeze(hrtf_meas(:,2,:));
    N = size(hrtf_l,2);
    [X,Y] = meshgrid( (1:N)',theta_measurement  );
    hrtf_interp = [ interp2(X,Y,hrtf_l,(1:N)',theta_interp*ones(N,1),'spline').';
         interp2(X,Y,hrtf_r,(1:N)',theta_interp*ones(N,1),'spline').'];

    case 'fourier'
%     % Interpolation in Fourier (wavenumber) domain
    hrtf_l = squeeze(hrtf_meas(:,1,:));
    hrtf_r = squeeze(hrtf_meas(:,2,:));
    L_f = fft(hrtf_l,[],1);
    R_f = fft(hrtf_r,[],1);
    K = size(theta_measurement,1);
    f = 1/K*exp(  1i*(0:K-1)'*theta_interp  );
    hrtf_interp = [(L_f.'*f).';(R_f.'*f).'];
end
hrtf_interp = squeeze(hrtf_interp);

end

