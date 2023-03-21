function setup = get_default_renderer_setup(renderer_type)
    switch renderer_type
        case 'Direct_playback'
            setup = [];
        case 'VBAP'
            setup = [];
        case 'DBAP'
            setup = [];
        case 'WFS'
            setup = struct('Antialiasing','on');
        case 'LWFS'
            setup = struct('Antialiasing','on');
        case 'NFC_HOA'
            setup = struct('M',30);
        case 'CTC'
            hrtf_sofa = SOFAload('BuK_ED_corr.sofa');
            setup = struct( 'Plant_model','HRTF','VS_model','HRTF', 'HRTF_database',hrtf_sofa,'N_filt',4096);
    end

end

