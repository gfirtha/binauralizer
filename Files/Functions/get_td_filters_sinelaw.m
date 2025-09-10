function [td_filters] = get_td_filters_sinelaw( xs, x0, n0, fs, c )
%xs = [2.5 0.5];
xr = [0 0];
N  = 512;
w = [ 0 : N/2 ]'/N*2*pi*fs;
v = bsxfun( @minus, xr, x0 );
v = bsxfun(@times, v,1./sqrt(sum(v.^2,2)));
vn = sum(v.*( (xr-xs)/norm(xr-xs) ),2);
[a,ix] = sort(vn,'descend');
R0 = mean(sqrt(sum(x0.^2,2)));

td_filters = zeros(N,size(x0,1));
if max(abs(vn)) == 1
    ixs = ix(1);
else
    ixs = ix(1:2);
    fi_base = acos((x0(ixs(1),:)-xr)/norm((x0(ixs(1),:)-xr))*...
        (x0(ixs(2),:)-xr)'/norm((x0(ixs(2),:)-xr)))/2;
    base_vec = sum(x0(ixs,:)-xr,1)/norm(sum(x0(ixs,:)-xr,1));
    base_vec = base_vec/norm(base_vec);
    fi0 = acos( base_vec*(xs-xr)'/norm((xs-xr)) );
%    H_sp1 = exp(1i* atan( w/c*R0*tan(fi0)/tan(fi_base) ) );
%    H_sp2 = exp(-1i*atan( w/c*R0*tan(fi0)/tan(fi_base) ) );
%    H_sp1(end) = real(H_sp1(end));
%    H_sp2(end) = real(H_sp2(end));
%     
    H_sp1_lin = exp(1i* ( 1/4*w/c*R0*(tan(fi0)/tan(fi_base) )) );
    H_sp2_lin = exp(-1i*( 1/4*w/c*R0*(tan(fi0)/tan(fi_base) )) );
    H_sp1_lin(end) = real(H_sp1_lin(end));
    H_sp2_lin(end) = real(H_sp2_lin(end));
%     
%      figure;
%      semilogx(w/(2*pi),20*log10(atan( w/c*R0*tan(fi0)/tan(fi_base) ))); 
%      hold on;
%      semilogx(w/(2*pi),20*log10( w/c*R0*tan(fi0)/tan(fi_base) ))
%      xlim([20,2e3])
%      
    td_filters(:,ixs(1)) = fftshift(ifft([H_sp1_lin;flipud(conj(H_sp1_lin(2:end-1)))]));
    td_filters(:,ixs(2)) = fftshift(ifft([H_sp2_lin;flipud(conj(H_sp2_lin(2:end-1)))]));
end

