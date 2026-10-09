classdef MasterEQ < handle
    % audioapp.dsp.MasterEQ
    % Single point of truth for master equalization.
    % - Enabled=false → no objects are allocated, apply() is a no-op
    % - Modes: 'graphic' (10-band) or 'lowpass'
    % - Toolboxes handled gracefully; graphicEQ needs Audio TB, lowpass needs DSP ST

    properties
        % config/state
        Enabled   (1,1) logical = false
        Mode      (1,:) char    = 'graphic'   % 'graphic' | 'lowpass'
        Fs        (1,1) double  = 48000

        % graphic EQ (Audio Toolbox)
        GraphicGains (1,10) double = zeros(1,10)
        EqGraphic     % graphicEQ or []

        % lowpass (DSP System Toolbox)
        LowpassFc (1,1) double = 20000
        EqLowpass   % dsp.LowpassFilter or []

        % environment flags
        HasDSP   (1,1) logical = ~isempty(ver('dsp'))
        HasAudio (1,1) logical = ~isempty(ver('audio'))
    end

    methods
        function obj = MasterEQ(ctxOrFs)
            if nargin>=1 && ~isempty(ctxOrFs)
                if isa(ctxOrFs,'audioapp.AppContext') || (isobject(ctxOrFs) && isprop(ctxOrFs,'Setup'))
                    try obj.Fs = ctxOrFs.Setup.SampleRate; catch, end
                    try obj.HasDSP   = ctxOrFs.HasDSP;   catch, end
                    try obj.HasAudio = ctxOrFs.HasAudio; catch, end
                elseif isnumeric(ctxOrFs)
                    obj.Fs = double(ctxOrFs);
                end
            end
        end

        function setFs(obj, fs)
            fs = max(1000, double(fs));
            if obj.Fs == fs, return; end
            obj.Fs = fs;
            if obj.Enabled
                obj.ensureObjects(); % rebuild with new Fs
            end
        end

        function setEnabled(obj, onOff)
            onOff = logical(onOff);
            if onOff == obj.Enabled, return; end
            obj.Enabled = onOff;
            if onOff
                obj.ensureObjects();
            else
                obj.destroyObjects();
            end
        end

        function setMode(obj, mode)
            mode = lower(char(mode));
            if ~ismember(mode, {'graphic','lowpass'}), return; end
            if strcmp(mode,'graphic') && ~obj.HasAudio
                % No Audio TB → force lowpass
                mode = 'lowpass';
            end
            obj.Mode = mode;
            if obj.Enabled
                obj.ensureObjects();
            end
        end

        function setGraphicGains(obj, gainsDB)
            if numel(gainsDB) ~= 10, return; end
            obj.GraphicGains = gainsDB(:).';
            if obj.Enabled && ~isempty(obj.EqGraphic) && isvalid(obj.EqGraphic)
                try, obj.EqGraphic.Gains = obj.GraphicGains; catch, obj.ensureObjects(); end
            end
        end

        function setLowpassFc(obj, fc)
            fc = max(50, min(0.499*obj.Fs, double(fc)));
            obj.LowpassFc = fc;
            if obj.Enabled && ~isempty(obj.EqLowpass)
                obj.ensureObjects(); % rebuild filter
            end
        end

        function setPreset(obj, name)
            name = lower(strtrim(string(name)));
            switch name
                case "flat"
                    obj.setGraphicGains(zeros(1,10));
                case "bass boost"
                    obj.setGraphicGains([+6 +4 +2 0 0 0 -1 -2 -3 -4]);
                case "treble cut"
                    obj.setGraphicGains([0 0 0 0 -2 -4 -6 -6 -6 -6]);
                case "vocal clarity"
                    obj.setGraphicGains([-2 -1 0 +2 +3 +3 +1 0 -1 -2]);
                case "lowpass 8k"
                    obj.setMode('lowpass'); obj.setLowpassFc(8000);
                case "lowpass 2k"
                    obj.setMode('lowpass'); obj.setLowpassFc(2000);
                otherwise
                    % no-op
            end
        end

        function y = apply(obj, x)
            % x: [Nsamp x Nchan]
            if ~obj.Enabled || isempty(x)
                y = x; return;
            end
            switch obj.Mode
                case 'graphic'
                    if ~isempty(obj.EqGraphic)
                        y = obj.EqGraphic(x);
                    else
                        y = x;
                    end
                case 'lowpass'
                    if ~isempty(obj.EqLowpass)
                        y = obj.EqLowpass(x);
                    else
                        y = x;
                    end
                otherwise
                    y = x;
            end
        end
    end

methods (Access=private)
    function ensureObjects(obj)
        obj.destroyObjects(); % clean slate

        switch obj.Mode
            case 'graphic'
                if obj.HasAudio && exist('graphicEQ','class') == 8
                    % Construct first, then set properties with fallbacks
                    G = obj.GraphicGains;
                    if numel(G) ~= 10, G = zeros(1,10); end
                    try
                        eq = graphicEQ();   % construct with defaults
                    catch ME
                        % If construction fails for any reason, fall back to low-pass
                        warning('MasterEQ:GraphicEQCreate', ...
                           'graphicEQ create failed (%s); falling back to low-pass.', ME.message);
                        obj.Mode = 'lowpass';
                        obj.ensureObjects();
                        return;
                    end

                    % Sample rate (if supported)
                    if isprop(eq,'SampleRate')
                        try, eq.SampleRate = obj.Fs; catch, end
                    end

                    % Number of bands varies by release: NumEQBands vs NumBands
                    nbProp = '';
                    if isprop(eq,'NumEQBands'), nbProp = 'NumEQBands'; end
                    if isempty(nbProp) && isprop(eq,'NumBands'), nbProp = 'NumBands'; end
                    if ~isempty(nbProp)
                        try
                            eq.(nbProp) = 10;
                        catch
                            % Some releases make it non-tunable; ignore if so
                        end
                    end

                    % Gains property name (usually 'Gains')
                    if isprop(eq,'Gains')
                        % Ensure length matches what object expects (query if available)
                        try
                            nReq = 10;
                            if ~isempty(nbProp)
                                nReq = eq.(nbProp);
                            end
                            if numel(G) ~= nReq
                                G = padarray(G(:).', [0 max(0,nReq-numel(G))], 0, 'post');
                                G = G(1:nReq);
                            end
                        catch
                        end
                        try, eq.Gains = G; catch, end
                    end

                    obj.EqGraphic = eq;
                else
                    % No Audio Toolbox → silently fall back
                    obj.Mode = 'lowpass';
                    obj.ensureObjects();
                end

            case 'lowpass'
                if obj.HasDSP
                    try, release(obj.EqLowpass); catch, end
                    % % Keep pass/stop close to requested Fc; clamp inside Nyquist
                    % Fc = max(50, min(0.499*obj.Fs, obj.LowpassFc));
                    % pb = max(40, 0.90*Fc);
                    % sb = max(45, Fc);
                    % obj.EqLowpass = dsp.LowpassFilter( ...
                    %     'SampleRate',          obj.Fs, ...
                    %     'FilterType',          'IIR', ...
                    %     'PassbandFrequency',   pb, ...
                    %     'StopbandFrequency',   sb, ...
                    %     'PassbandRipple',      0.5, ...
                    %     'StopbandAttenuation', 60);

                    Fc = max(50, min(0.499*obj.Fs, obj.LowpassFc));

                    TW = max( max(60, 0.6*Fc),  0.10*obj.Fs );
                    pb = max(40, Fc - TW/2);
                    sb = min(0.499*obj.Fs, Fc + TW/2);

                    Rp = 0.1;
                    As = 60; 

                    d = designfilt('lowpassfir', ...
                        'PassbandFrequency',   pb, ...
                        'StopbandFrequency',   sb, ...
                        'PassbandRipple',      Rp, ...
                        'StopbandAttenuation', As, ...
                        'DesignMethod',        'kaiserwin', ...
                        'SampleRate',          obj.Fs);
                    obj.EqLowpass = dsp.FIRFilter('Numerator', d.Coefficients);

                else
                    obj.EqLowpass = [];
                end
        end
    end

    function destroyObjects(obj)
        try, release(obj.EqLowpass); catch, end
        obj.EqLowpass = [];
        obj.EqGraphic = [];
    end
end

end
