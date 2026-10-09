function directivity_mx = get_twoway_dir(fs,N_fft,fi,R)
w = [(0 : N_fft/2 - 1)';(-N_fft/2:-1)' ]/N_fft*2*pi*fs;
c = audioapp.util.Physics.speedOfSound();

fc = 400;
wc = fc*2*pi;
[Fi, W] = meshgrid(fi,w);

Nfilter = 4;
LP = (1./sqrt(1+(W/wc).^(2*Nfilter))).^2;
HP = (1./sqrt(1+(wc./W).^(2*Nfilter))).^2;

D = 2*(besselj(1,R(1)/c*W.*sin(Fi))./(R(1)/c*W.*sin(Fi)).*LP*1+...
       besselj(1,R(2)/c*W.*sin(Fi))./(R(2)/c*W.*sin(Fi)).*HP*1).*exp(-1i*W*Nfilter/fs/2);
D(1,:) = 1;
D(:,1) = 1;
D(end/2+1,:) = real(D(end/2+1,:));
directivity_mx = D;
end

