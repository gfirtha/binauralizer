function [position, normal] = get_default_layout(Num_channels, R, Shape)

switch Shape
    case 'circular'
        switch Num_channels
            case 1
                [x,y] = pol2cart(0*pi/180 + pi/2,R);
            case 2
                [x,y] = pol2cart([30; -30]*pi/180 + pi/2,repmat(R,Num_channels,1));
            case 3
                [x,y] = pol2cart([30; -30; 0]*pi/180 + pi/2,repmat(R,Num_channels,1));
            case 4
                %            [x,y] = pol2cart([45; -45; -135; -225]*pi/180,repmat(R,Num_channels,1));
                [x,y] = pol2cart([-30; 0; 30; 180]*pi/180 + pi/2,repmat(R,Num_channels,1));
            case 5
                [x,y] = pol2cart([30; 0; -30; 110; -110]*pi/180 + pi/2,repmat(R,Num_channels,1));
            case 6
                [x,y] = pol2cart([-30; 0; 30; 110; -110]*pi/180 + pi/2,repmat(R,Num_channels-1,1));
                x = [x;0.2];
                y = [y;0.2];
            otherwise
                [x,y] = pol2cart(linspace(0,360*(1-1/Num_channels), Num_channels)'*pi/180 + pi/2,repmat(R,Num_channels,1));
        end
        position = [x,y];
        normal = -position./sqrt(sum(position.^2,2));
    case 'linear'
        position = [linspace(-R,R,Num_channels)',R*ones(Num_channels,1)];
        normal = [ zeros(Num_channels,1), -ones(Num_channels,1) ];
    case 'stereo'
        [x,y] = pol2cart([30; -30]*pi/180 + pi/2,repmat(R,2,1));
        position = [x,y];
        normal = -position./sqrt(sum(position.^2,2));
    case '5.0'
        [x,y] = pol2cart([30; 0; -30; 110; -110]*pi/180 + pi/2,repmat(R,5,1));
        position = [x,y];
        normal = -position./sqrt(sum(position.^2,2));
    case '5.1'
        [x,y] = pol2cart([-30; 0; 30; 110; -110]*pi/180 + pi/2,repmat(R,Num_channels-1,1));
        x = [x;0.2];
        y = [y;0.2];
        position = [x,y];
        normal = -position./sqrt(sum(position.^2,2));
    case 'general'

        % default folder "SSDs"
        startDir = fullfile(pwd,'Data','SSDs');
        if ~exist(startDir,'dir'), startDir = pwd; end
        [f,p] = uigetfile({'*.eps','EPS files (*.eps)'}, 'Select SSD geometry (EPS)', startDir);
        if isequal(f,0), return; end
        fn = fullfile(p,f);

        % N comes from UI / setup
        N = Num_channels;
        if isnan(N) || N<=2 || mod(N,1)~=0
            warndlg('Please set a valid integer LS number (N > 2) first.','Invalid N','modal');
            return;
        end

        % Parse EPS → sample N points → inward normals
        try
            [x0,n0] = secondary_source_distribution.load_eps_as_ssd(fn, N,R);
        catch ME
            errordlg(sprintf('Failed to load geometry:\n\n%s', ME.message),'Load failed','modal');
            return;
        end
        
        position = x0;
        normal = n0;



end

end

