function [mx] = get_downmixing_mx(loudspeakers,Nout)

if (length(loudspeakers) <= Nout)
    mx = eye(length(loudspeakers));
else
    ls_in = cell2mat(cellfun( @(x) x.position', loudspeakers,'UniformOutput',false));
    R0 = cell2mat(cellfun( @(x) sqrt(sum(x.position.^2,2)), loudspeakers,'UniformOutput',false));
    ls_in = bsxfun(@times, ls_in, 1./R0);
    fi_in = cart2pol(ls_in(1,:),ls_in(2,:));
    ls_out = get_default_layout(Nout, 1)';
    fi_out = cart2pol(ls_out(1,:),ls_out(2,:));

    for n = 1 : length(fi_in)
        [~,i] = sort(abs( fi_out - fi_in(n) ));
        l = ls_out(:,sort(i(1:2)));
        ind = i(1:2);
        g = l\ls_in(:,n)/norm(l\ls_in(:,n));
        mx(:,n) = g/norm(g);
    end
    

end
end

