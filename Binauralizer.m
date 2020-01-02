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
SOFAstart;
hrtf_sofa = SOFAload('FABIAN_HRIR_measured_HATO_0.sofa');
%hrtf_sofa = SOFAload('subject_003.sofa');
%hrtf_sofa = SOFAload('C1m.sofa');

handles.binauralizer_setup = struct(   ...
    'Input_file',               'eli_50.wav',...
    'Block_size',               1024*2, ...
    'HRTF',                     hrtf_sofa, ...
    'Volume',                   0.5, ...
    'Rendering',                'Binaural',...
    'WFS_Setup',                struct('R',2,'N',32));

handles.fileReader = dsp.AudioFileReader(handles.binauralizer_setup.Input_file,...
                       'SamplesPerFrame',handles.binauralizer_setup.Block_size,'PlayCount',10);
handles.gui = listener_space_axes;
handles.gui.draw_gui(handles.axes1);
handles.sound_scene = sound_scene(handles.fileReader,handles.gui,handles.binauralizer_setup);
handles.output = hObject; 
guidata(hObject, handles);
handles

% --- Outputs from this function are returned to the command line.
function varargout = Binauralizer_OutputFcn(hObject, eventdata, handles)
warning('off','all')
varargout{1} = handles.output;

% --- Executes on button press in play_btn.
function play_btn_Callback(hObject, eventdata, handles)
deviceWriter = audioDeviceWriter('SampleRate',handles.fileReader.SampleRate);
handles.stop_now = 0;
guidata(hObject,handles);
elapsed_time = 0;
while (~isDone(handles.fileReader))&&(~handles.stop_now)
    output = handles.sound_scene.binauralize_sound_scene(handles.binauralizer_setup.Volume*handles.fileReader());
    deviceWriter(output);
    drawnow limitrate
    handles = guidata(hObject);
end
release(deviceWriter)
release(handles.fileReader)
fprintf('\n Elapsed time: %d \n',mean(elapsed_time))

% --- Executes on button press in load_file_btn.
function load_file_btn_Callback(hObject, eventdata, handles)
[file,path] = uigetfile('*.wav;*.mp3;*.aac;*.ac3');
handles.binauralizer_setup.Input_file = strcat(path,file);
handles.fileReader = dsp.AudioFileReader(handles.binauralizer_setup.Input_file,...
                       'SamplesPerFrame',handles.binauralizer_setup.Block_size);
audio_info = info(handles.fileReader);
handles.binauralizer_setup.Default_source_position = get_default_layout(audio_info.NumChannels,1.5);
for n = 1 : length(handles.source_points)
    delete(handles.source_points{n});
end
[handles.source_points, handles.source_listeners] = handles.gui.add_sources(handles.binauralizer_setup,handles.axes1);
guidata(hObject,handles);

% --- Executes on slider movement.
function Volume_Callback(hObject, eventdata, handles)
handles.binauralizer_setup.Volume = get(hObject,'Value');
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Volume_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in stop_btn.
function stop_btn_Callback(hObject, eventdata, handles)
handles.stop_now = 1;
guidata(hObject, handles);

% TODO: Warning due to drawpoint
%function allevents(src,evt)
%evname = evt.EventName;