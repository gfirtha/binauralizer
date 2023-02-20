function directivity_mx = get_piston_dir(fs,N_fft,fi,R)
w = [(0 : N_fft/2 - 1)';(-N_fft/2:-1)' ]/N_fft*2*pi*fs;
c = 343.1;
k = w / c;
[Fi, K] = meshgrid(fi,k);
D = 2*besselj(1,R*K.*sin(Fi))./(R*K.*sin(Fi));
D(1,:) = 1;
D(:,1) = 1;
D(end/2+1,:) = real(D(end/2+1,:));
directivity_mx = D;
end

