classdef VideoService < handle
    properties
        C
        Reader
        FPS double = NaN
        SpfProc double = NaN
        Accum double = 0
        Fig; Ax; Im
        Title char = 'Video preview'
    end

    methods
        function obj = VideoService(ctx), obj.C = ctx; end
        function tf = hasVideo(obj), tf = ~isempty(obj.Reader); end

        function open(obj, filename)
            obj.close();
            obj.Reader = VideoReader(filename);
            obj.FPS    = obj.Reader.FrameRate;
            obj.SpfProc= obj.C.Setup.SampleRate / max(1e-9, obj.FPS);
            [~,bn,ext] = fileparts(filename);
            obj.Title  = sprintf('Video: %s%s', bn, ext);
            obj.Accum  = 0;
        end

        function rewind(obj)
            if obj.hasVideo()
                try, obj.Reader.CurrentTime = 0; catch, end
                obj.Accum = 0;
            end
        end

        function ensureWindow(obj)
            if isempty(obj.Fig) || ~ishandle(obj.Fig)
                obj.Fig = figure('Name',obj.Title,'NumberTitle','off','MenuBar','none','ToolBar','none',...
                                 'Color',[0 0 0],'Units','pixels','Position',[100 100 960 540],...
                                 'CloseRequestFcn',@(~,~)obj.onClosed());
                obj.Ax = axes('Parent',obj.Fig); axis(obj.Ax,'image','off');
                obj.Im = [];
            end
        end

        function showFirstFrame(obj)
            if ~obj.hasVideo(), return; end
            obj.ensureWindow();
            try, obj.Reader.CurrentTime = 0; catch, end
            obj.Accum = 0;
            if hasFrame(obj.Reader)
                f0 = readFrame(obj.Reader);
                if isempty(obj.Im) || ~ishandle(obj.Im)
                    obj.Im = imshow(f0,'Parent',obj.Ax);
                else
                    set(obj.Im,'CData',f0);
                end
                drawnow;
            end
        end

        function onAudioSamplesProcessed(obj, nSamples)
            if ~obj.hasVideo() || isempty(obj.Im) || ~ishandle(obj.Im)
                return;
            end
            obj.Accum = obj.Accum + nSamples;
            while (obj.Accum >= obj.SpfProc) && hasFrame(obj.Reader)
                obj.Accum = obj.Accum - obj.SpfProc;
                f = readFrame(obj.Reader);
                set(obj.Im,'CData',f);
            end
            drawnow limitrate
        end

        function onClosed(obj, varargin)
            if ~isempty(obj.Fig) && ishandle(obj.Fig), delete(obj.Fig); end
            obj.Fig = []; obj.Ax = []; obj.Im = [];
        end

        function close(obj)
            obj.onClosed();
            obj.Reader = [];
            obj.FPS = NaN; obj.SpfProc = NaN; obj.Accum = 0;
        end
    end
end
