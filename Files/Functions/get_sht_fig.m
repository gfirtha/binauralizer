function [N,M,out] = get_sht_fig(Nmax,fix,ch, input)

[N,M] = meshgrid((0:Nmax),(-Nmax:Nmax));
nlin = N.^2+N+M+1;
nlin(abs(M)>N) = 0;

out = interp1( (1:size(input,1)), squeeze(input(:,ch,fix)), nlin,'linear');

end

