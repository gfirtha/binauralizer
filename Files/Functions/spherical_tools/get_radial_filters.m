function [extrapolation_filters] = get_radial_filters(r,R0,Nf,fs,Nmax)
type = 'PW_approx';
c = 343.1;
omega = 2*pi*(0:Nf-1)'/Nf*fs;

switch type
    case 'physical'
        %Physical solution
        Hr  = getSphH((0:Nmax)',2,omega/c*r);
        Hr0 = getSphH((0:Nmax)',2,omega/c*R0);
        H_radial = Hr./Hr0;
        if r < R0
            fcs = (1:Nmax)./R0*c/pi/1.2;
            Hreg = get_regularization_filters( omega, fcs );
            H_radial = H_radial.*Hreg.';
           %sd figure;semilogx(omega/2/pi,20*log10(abs(H_radial.')));ylim([-20,10]);grid on
        end
    case 'PW_approx'
        % High-freq. approx
        Hr  = (-1i).^(0:Nmax)'.*(exp(-1i*omega/c*r)./(omega/c*r)).';
        Hr0  = (-1i).^(0:Nmax)'.*(exp(-1i*omega/c*R0)./(omega/c*R0)).';
        H_radial = Hr./Hr0;
end
H_radial(isnan(H_radial)) = 0;

for n = 0 : Nmax
    for m = -n:n
        extrapolation_filters(n^2+n+m+1,:) = H_radial(n+1,:);
    end
end

end




