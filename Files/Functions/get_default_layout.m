function [position, normal, channelMap] = get_default_layout(Num_channels, R, Shape, opts)
%GET_DEFAULT_LAYOUT  Loudspeaker positions, normals, and channel mapping.
%
%   [position, normal] = get_default_layout(N, R, Shape)
%   [position, normal, channelMap] = get_default_layout(N, R, Shape)
%   [...] = get_default_layout(N, R, Shape, opts)   opts: Loudspeaker_setup
%                (layout options, e.g. Surround51 for 'soundbar')
%
%   position   : N x 2, loudspeaker coordinates (SSD contour order)
%   normal     : N x 2, inward-facing unit normals
%   channelMap : 1 x N, SSD index → physical output channel index.
%                For standard layouts (5.0, 5.1, 7.1), this maps the
%                geometrical (contour) order to the ITU/SMPTE channel order.
%                For linear/circular arrays, channelMap = 1:N (identity).
%
%   The SSD contour order is chosen so that loudspeakers are sequential
%   along the contour (e.g., counter-clockwise). This is what all renderers
%   (WFS, HOA, VBAP, ...) expect internally.
%
%   The channel order is the physical output channel assignment:
%     5.0 ITU:  ch1=L(30°)  ch2=R(-30°)  ch3=C(0°)  ch4=LS(110°)  ch5=RS(-110°)
%     5.1 ITU:  ch1=L  ch2=R  ch3=C  ch4=LFE  ch5=LS  ch6=RS

channelMap = 1:Num_channels;   % default: identity

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
                [x,y] = pol2cart([-30; 0; 30; 180]*pi/180 + pi/2,repmat(R,Num_channels,1));
            case 5
                % SSD contour order (CCW from front-left):
                %   SSD#1 = L(30°), SSD#2 = C(0°), SSD#3 = R(-30°),
                %   SSD#4 = LS(110°), SSD#5 = RS(-110°)
                [x,y] = pol2cart([30; 0; -30; 110; -110]*pi/180 + pi/2,repmat(R,Num_channels,1));
                % ITU 5.0: ch1=L, ch2=R, ch3=C, ch4=LS, ch5=RS
                % SSD#1→ch1(L), SSD#2→ch3(C), SSD#3→ch2(R), SSD#4→ch4(LS), SSD#5→ch5(RS)
                channelMap = [1, 3, 2, 4, 5];
            case 6
                % SSD order: L, C, R, LS, RS, LFE.
                [x,y] = pol2cart([30; 0; -30; 110; -110]*pi/180 + pi/2,repmat(R,Num_channels-1,1));
                x = [x;0.2];
                y = [y;0.2];
                channelMap = [1, 3, 2, 5, 6, 4];
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
        % SSD contour order (front-left, center, front-right, left-surr, right-surr)
        [x,y] = pol2cart([30; 0; -30; 110; -110]*pi/180 + pi/2,repmat(R,5,1));
        position = [x,y];
        normal = -position./sqrt(sum(position.^2,2));
        % ITU 5.0: ch1=L, ch2=R, ch3=C, ch4=LS, ch5=RS
        channelMap = [1, 3, 2, 4, 5];

    case '5.1'
        [x,y] = pol2cart([30; 0; -30; 110; -110]*pi/180 + pi/2,repmat(R,Num_channels-1,1));
        x = [x;0.2];
        y = [y;0.2];
        position = [x,y];
        normal = -position./sqrt(sum(position.^2,2));
        % ITU 5.1: ch1=L, ch2=R, ch3=C, ch4=LFE, ch5=LS, ch6=RS
        % SSD: #1=L, #2=C, #3=R, #4=LS, #5=RS, #6=LFE
        channelMap = [1, 3, 2, 5, 6, 4];

    case {'5.1.4', '7.1.4'}
        % Dolby Atmos home layouts (Dolby Atmos Home Theater Installation
        % Guidelines): ear-level speakers at z = 0, the four overhead
        % speakers at +-45 deg (top front) / +-135 deg (top rear) azimuth,
        % 45 deg elevation, all on a sphere of radius R around the listener.
        % No loudspeaker for the LFE: it has no direction; its output
        % channel (4) stays reserved, so the channel order is the
        % WAVE/Windows one (WAVEFORMATEXTENSIBLE speaker mask order):
        %   5.1.4: FL FR FC LFE SL SR TFL TFR TBL TBR            (10 ch)
        %   7.1.4: FL FR FC LFE BL BR SL SR TFL TFR TBL TBR      (12 ch)
        % SSD order: ear level counter-clockwise from front-left, then the
        % height layer counter-clockwise.  [azimuth (CCW, 0 = front), elevation]
        if strcmp(Shape, '5.1.4')
            dirs = [30 0; 0 0; -30 0; -110 0; 110 0; ...       % L C R  Rs Ls
                    45 45; -45 45; -135 45; 135 45];           % TFL TFR TBR TBL
            channelMap = [1 3 2 6 5, 7 8 10 9];
        else
            dirs = [30 0; 0 0; -30 0; -90 0; -150 0; 150 0; 90 0; ...  % L C R Rss Rrs Lrs Lss
                    45 45; -45 45; -135 45; 135 45];
            channelMap = [1 3 2 8 6 5 7, 9 10 12 11];
        end
        az = dirs(:,1)*pi/180 + pi/2;          % +y is the front
        el = dirs(:,2)*pi/180;
        position = R*[cos(el).*cos(az), cos(el).*sin(az), sin(el)];
        normal = -position./sqrt(sum(position.^2,2));
        Num_channels = size(dirs,1);

    case 'soundbar'
        % The 16-driver WFS soundbar (2 x 8-driver modules, 62.5 mm pitch,
        % 1 m long), ear height, exactly between the front speakers: at
        % depth R*cos(30 deg) (1.3 m for the default R = 1.5 m).
        % Optional discrete surrounds (opts.Surround51): a standard 5.1
        % without the centre - the bar is the centre - fronts at +-30 deg
        % on radius R, surrounds at +-110 deg on the closer radius 0.8 R.
        % Output channels: bar 1..16, then L R LFE Ls Rs (17..21, LFE
        % reserved, no loudspeaker).
        Nbar = 16; pitch = 0.0625;
        x = ((1:Nbar)' - (Nbar+1)/2)*pitch;
        position = [x, R*cosd(30)*ones(Nbar,1)];
        normal = [zeros(Nbar,1), -ones(Nbar,1)];
        channelMap = 1:Nbar;
        if nargin > 3 && isfield(opts, 'Surround51') && opts.Surround51
            az = [30; -30; 110; -110]*pi/180 + pi/2;      % L R Ls Rs
            [xs, ys] = pol2cart(az, R*[1; 1; 0.8; 0.8]);
            position = [position; xs, ys];
            normal = [normal; -[xs, ys]./sqrt(xs.^2 + ys.^2)];
            channelMap = [channelMap, Nbar + [1 2 4 5]];
        end
        Num_channels = size(position,1);

    case 'genius_array'
        % 8-element linear array, 7 cm inter-element spacing
        % Centered at y = R (like a standard linear array)
        Nel = 8;
        d_elem = 0.07;   % 7 cm
        total_width = (Nel-1) * d_elem;
        x = linspace(-total_width/2, total_width/2, Nel)';
        y = R * ones(Nel, 1);
        position = [x, y];
        normal = [zeros(Nel,1), -ones(Nel,1)];   % pointing inward (towards y=0)
        channelMap  = [8, 7, 4, 3, 6, 5, 2, 1];
        Num_channels = Nel;  % override — genius_array always has 8 elements

    case 'general'
        % default folder "SSDs"
        startDir = fullfile(pwd,'Data','SSDs');
        if ~exist(startDir,'dir'), startDir = pwd; end
        [f,p] = uigetfile({'*.eps','EPS files (*.eps)'}, 'Select SSD geometry (EPS)', startDir);
        if isequal(f,0), return; end
        fn = fullfile(p,f);

        N = Num_channels;
        if isnan(N) || N<=2 || mod(N,1)~=0
            warndlg('Please set a valid integer LS number (N > 2) first.','Invalid N','modal');
            return;
        end

        try
            [x0,n0] = secondary_source_distribution.load_eps_as_ssd(fn, N,R);
        catch ME
            errordlg(sprintf('Failed to load geometry:\n\n%s', ME.message),'Load failed','modal');
            return;
        end
        
        position = x0;
        normal = n0;
        channelMap = 1:N;   % custom geometries: identity by default
end

% ── Ensure 3D: pad position and normal to N×3 if needed ──
if size(position, 2) < 3
    position = [position, zeros(size(position,1), 3-size(position,2))];
end
if size(normal, 2) < 3
    normal = [normal, zeros(size(normal,1), 3-size(normal,2))];
end

end
