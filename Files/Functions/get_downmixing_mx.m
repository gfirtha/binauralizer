function [mx] = get_downmixing_mx(loudspeakers,Nout)

if (length(loudspeakers) <= Nout)
    mx = eye(length(loudspeakers));
else
    ls_in = cell2mat(cellfun( @(x) x.position', loudspeakers,'UniformOutput',false));
    R0 = cell2mat(cellfun( @(x) sqrt(sum(x.position.^2,2)), loudspeakers,'UniformOutput',false));
    ls_in = bsxfun(@times, ls_in, 1./R0);
    
    Nout = min(Nout,5);
    ls_out = get_default_layout(Nout, 1, 'circular')';

    mx = zeros(size(ls_out,2),size(ls_in,2));
    for n = 1 : size(ls_in,2)
        [~,i] = sort(ls_out'*ls_in(:,n),'descend');
        ind = i(1:2);
        l = ls_out(:,ind);
        g = l\ls_in(:,n)/norm(l\ls_in(:,n));
        mx(ind,n) = g/norm(g);
    end
    

end
end

