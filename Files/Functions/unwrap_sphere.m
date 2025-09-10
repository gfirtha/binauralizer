function [phase] = unwrap_sphere(sphere,phase,ix_start)
N0 = 100;
figure
pl = sphere.plot_data(squeeze(phase(:,N0)));
axis equal tight

adjacency = sphere.get_delanuey_adjacency;

R = get_rotation_mx(sphere.x0(ix_start,:)',[0,0,1]');
Xout = (R*sphere.x0')';
[~,elev] = cart2sph(Xout(:,1),Xout(:,2),Xout(:,3));
zenith = pi/2-elev;

ixs_to_unwrap = setdiff((1:size(sphere.x0,1)),ix_start);
ixs_unwrapped = ix_start;

while ~isempty(ixs_to_unwrap)

    neighbours_to_unwrap = setdiff(adjacency{ix_start},ixs_unwrapped);
    if isempty(neighbours_to_unwrap)
        [~,ixmin] = min(sum((sphere.x0(ix_start,:)-sphere.x0(ixs_to_unwrap,:)).^2,2));
        ix_curr = ixs_to_unwrap(ixmin);
    else
        [~,ixmin] = min(zenith(neighbours_to_unwrap));
        ix_curr = neighbours_to_unwrap(ixmin);
    end

    ixs = intersect(adjacency{ix_curr},ixs_unwrapped);
%    phase(ix_curr,:) = min(phase(ixs,:) + mod(phase(ix_curr,:)-phase(ixs,:)+pi,2*pi)-pi,[],1);
    
   [~,ixmin] = min(abs(mod(phase(ix_curr,:)-phase(ixs,:)+pi,2*pi)-pi),[],1);
   ix0 = sub2ind(size(phase),ixs(ixmin)', (1:size(phase,2))');
   phase(ix_curr,:) = phase(ix0)' + mod(phase(ix_curr,:)-phase(ix0)'+pi,2*pi)-pi;
   

    ixs_unwrapped = union(ixs_unwrapped, ix_curr);
    ixs_to_unwrap = setdiff(ixs_to_unwrap,ix_curr);
    set(pl,'CData',phase(:,100))
    drawnow

    ix_start = ix_curr;
end

end