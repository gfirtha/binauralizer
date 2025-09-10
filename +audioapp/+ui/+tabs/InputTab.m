classdef InputTab < handle
    properties
        C
        Panel
        BtnLoad    matlab.ui.control.UIControl
        Table      audioapp.ui.SourceTable
    end

    methods
        function obj = InputTab(parentTab, ctx)
            obj.C = ctx;
            obj.Panel = parentTab;

            obj.BtnLoad = uicontrol(obj.Panel,'Style','pushbutton','String','Load media ...',...
                'Units','normalized','Position',[0.07 0.88 0.86 0.08],...
                'Callback',@(s,e)obj.onLoadMedia());

            obj.Table = audioapp.ui.SourceTable(obj.Panel, [0.02 0.06 0.975 0.78], ctx);
        end

        function onLoadMedia(obj)
            startDir = fullfile(pwd,'Data','SoundSamples');
            [f,p] = uigetfile( ...
                {'*.wav;*.mp3;*.aac;*.ac3;*.flac;*.m4a;*.ogg;*.wma;*.mp4;*.mov;*.m4v;*.avi', ...
                 'Audio or video files'}, ...
                'Select media',startDir);
            if isequal(f,0), return; end
            obj.C.ensureEngine();
            kind = obj.C.Engine.loadMedia(fullfile(p,f)); %#ok<NASGU>
            obj.C.rebuildSoundScene();
            if obj.C.Video.hasVideo(), obj.C.Video.showFirstFrame(); end
            obj.Table.refresh();
        end
    end
end
