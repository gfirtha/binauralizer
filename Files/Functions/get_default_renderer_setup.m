function setup = get_default_renderer_setup(renderer_type)
    switch renderer_type
        case 'Direct_playback'
            setup = [];
        case 'VBAP'
            setup.VBAPspread  = 0;    % spread angle (degrees)
            setup.VBAPpNorm   = 2;    % p-norm: 1 = amplitude, 2 = power
            setup.VBAPnSpread = 8;    % number of auxiliary spread sources
        case 'DBAP'
            setup = [];
        case 'TD_stereo'
            setup.TDNumTaps     = 2048;
            setup.TDMethod      = 'firtha_linear';
            setup.TDRefPoint    = [0 0 0];
            setup.TDRefRadiusR0 = 'auto';
            setup.TDHybridAlpha = 0.5;
        case 'WFS'
            setup = struct('Antialiasing','off');
        case 'LWFS'
            setup = struct('Antialiasing','on');
        case 'LWFS_foc'
            setup = struct('Antialiasing','on');
        case 'LWFS_nara'
            setup = struct('Antialiasing','on');
        case 'Nara'
            setup = [];
        case 'HOA'
            setup = struct('M',20);
        case 'NFC_HOA'
            setup = struct('M',50);
        case 'CTC'
            hrtf_sofa = SOFAload('BuK_ED_corr.sofa');
    %       setup = struct( 'Plant_model','point_source','VS_model','point_source', 'HRTF_database',hrtf_sofa,'N_filt',1024);
     %       setup = struct( 'Plant_model','rigid_sphere','VS_model','rigid_sphere', 'HRTF_database',hrtf_sofa,'N_filt',1024);
            setup = struct( 'Plant_model','HRTF','VS_model','HRTF', 'HRTF_database',hrtf_sofa,'N_filt',1024);
        case 'Dolby_surround'
            setup = [];
        case 'Room_simulator'
    end

end

