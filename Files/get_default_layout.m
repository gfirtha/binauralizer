function [position] = get_default_layout(Num_channels, R)

    switch Num_channels
        case 1
            [x,y] = pol2cart(0*pi/180,R);
        case 2
            [x,y] = pol2cart([30; -30]*pi/180,repmat(R,Num_channels,1));
        case 3
            [x,y] = pol2cart([30; -30; 0]*pi/180,repmat(R,Num_channels,1));
        case 5
            [x,y] = pol2cart([30; -30; 0; 110; -110]*pi/180,repmat(R,Num_channels,1));
        otherwise
            [x,y] = pol2cart(linspace(0,360, Num_channels)'*pi/180,repmat(R,Num_channels,1));
    end
    position = [x,y];

end

