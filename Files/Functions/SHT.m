function [ out ] = SHT( varargin )
% A numerikus SHT megvalósítása
%in, zenith, azim, N, mode, type

in = varargin{1};
if length(varargin) == 6
    zenith = varargin{2};
    azim = varargin{3};
    N = varargin{4};
    mode = varargin{5};
    type = varargin{6};
else
    sphere = varargin{2};
    N = varargin{3};
    mode = varargin{4};
    type = varargin{5};
end

switch mode
    case 'trans'
        Y = getSpherHarmMx( zenith, azim, N, type );
        out = zeros(size(Y,2),size(in,2),size(in,3));
        for m = 1 : size(in,3)
            out(:,:,m) = (4*pi/length(zenith)*Y')*squeeze(in(:,:,m));
        end
    case 'pinv'
        Y = getSpherHarmMx( zenith, azim, N, type );
        out = zeros(size(Y,2),size(in,2),size(in,3));
        Yi = (Y'*Y)\Y';
        for m = 1 : size(in,3)
            out(:,:,m) = Yi*squeeze(in(:,:,m));
        end
    case 'num_int'
        quad_in = sphere.interpolate2quadrature(in);
        [azim, elev]= cart2sph(sphere.quadrature_points.quad_pos(:,1),sphere.quadrature_points.quad_pos(:,2),sphere.quadrature_points.quad_pos(:,3));
        zenith = pi/2 - elev;
        S = 4*pi*mean(sphere.Phi0(:,3))^2;
        Y = getSpherHarmMx( zenith, azim, N, type );
        out = zeros(size(Y,2),size(in,2),size(in,3));
        for m = 1 : size(in,3)
            out(:,:,m) = (4*pi*sphere.quadrature_points.quad_w.*Y./S)'*squeeze(quad_in(:,:,m));
        end

end

end

