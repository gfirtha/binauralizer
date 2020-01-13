function varargout = Binauralizer(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Binauralizer_OpeningFcn, ...
    'gui_OutputFcn',  @Binauralizer_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% --- Executes just before Binauralizer is made visible.
function Binauralizer_OpeningFcn(hObject, eventdata, handles, varargin)
addpath('Samples')
addpath(genpath('Files'))
addpath(genpath('RIRs'))
SOFAstart;
hrtf_sofa = SOFAload('FABIAN_HRIR_measured_HATO_0.sofa');
%hrtf_sofa = SOFAload('subject_003.sofa');
%hrtf_sofa = SOFAload('C1m.sofa');

% Setup options:
%   Rendering:              'Binaural':        Binaural rendering only 
%                           'WFS':             Wave Field Synthesis driving
%                                              functions 
%                           'VBAP':            Vector-based Amplitude
%                                              Panning
%   Binaural_source_type:   'point_source':    omnidirectional directivity,
%                                              parameter: R (only for drawing)
%                           'circular_piston': baffled piston directivity, 
%                                              parameter: R radius 
%                           'two_way_speaker': two baffled pistons, with
%                                              Linkwitz-Riley crossover filter,
%                                              parameter: [R_lp, R_hp] radii
%
%
%
%
handles.sound_scene_setup = struct(  ...
    'Input_file',               'gitL.wav',...
    'Block_size',               1024*4, ...
    'HRTF',                     hrtf_sofa, ...
    'Volume',                   0.5, ...
    'Binaural_source_type',     struct('Shape','circular_piston','R',0.06),...
    'Virtual_source_type',      struct('Shape','point_source','R',0.05),...
    'Rendering',                'WFS',...
    'renderer_setup',           struct('R',2,'N',128/2,'Antialiasing','on'));

handles.sound_scene_setup.Input_stream = dsp.AudioFileReader(handles.sound_scene_setup.Input_file,...
                       'SamplesPerFrame',handles.sound_scene_setup.Block_size,'PlayCount',10);
handles.gui = listener_space_axes(handles.axes1);
handles.sound_scene = sound_scene(handles.gui,handles.sound_scene_setup);
handles.output = hObject; 
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = Binauralizer_OutputFcn(hObject, eventdata, handles)
warning('off','all')
varargout{1} = handles.output;

% --- Executes on button press in play_btn.
function play_btn_Callback(hObject, eventdata, handles)
handles.stop_now = 0;
deviceWriter = audioDeviceWriter('SampleRate',handles.sound_scene_setup.Input_stream.SampleRate);
guidata(hObject,handles);
elapsed_time = 0;
i = 1;
while (~isDone(handles.sound_scene_setup.Input_stream))&&(~handles.stop_now)
    tic
    output = handles.sound_scene.binauralize_sound_scene(handles.sound_scene_setup.Volume*...
                                                         handles.sound_scene_setup.Input_stream());
    elapsed_time(i) = toc;
    deviceWriter(output);
    drawnow limitrate
    handles = guidata(hObject);
    i = i + 1;
   % mean(elapsed_time)
end
release(deviceWriter)
release(handles.sound_scene_setup.Input_stream)

% --- Executes on button press in load_file_btn.
function load_file_btn_Callback(hObject, eventdata, handles)
[file,path] = uigetfile('*.wav;*.mp3;*.aac;*.ac3');
handles.sound_scene_setup.Input_file = strcat(path,file);
handles.sound_scene_setup.Input_stream = dsp.AudioFileReader(handles.sound_scene_setup.Input_file,...
                       'SamplesPerFrame',handles.sound_scene_setup.Block_size);
handles.sound_scene.delete(handles.gui);
handles.sound_scene = sound_scene(handles.gui,handles.sound_scene_setup);
guidata(hObject,handles);

% --- Executes on slider movement.
function Volume_Callback(hObject, eventdata, handles)
handles.sound_scene_setup.Volume = get(hObject,'Value');
guidata(hObject,handles)

function Volume_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function stop_btn_Callback(hObject, eventdata, handles)
handles.stop_now = 1;
guidata(hObject, handles);