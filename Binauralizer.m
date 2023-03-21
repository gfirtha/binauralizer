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
install;
hrtf_sofa = SOFAload('BuK_ED_corr.sofa');
input_file = 'GitL.wav';
block_size = 1024;
handles.Volume = 0.5;

audiodevreset;
handles.Input_stream = dsp.AudioFileReader(input_file,...
    'SamplesPerFrame',block_size,'PlayCount',10);
handles.Output_device = audioDeviceWriter;
set(handles.Output_device,'SampleRate',handles.Input_stream.SampleRate);
set(handles.Output_device,'BufferSize',block_size);

Nch_in = handles.Input_stream.info.NumChannels;
Nch_out = handles.Output_device.info.MaximumOutputChannels;
fs = handles.Input_stream.SampleRate;

handles.sound_scene_setup = struct(  ...
    'N_in',                     Nch_in,...
    'N_out',                    Nch_out,...
    'SampleRate',               fs,...
    'Block_size',               block_size, ...
    'HRTF',                     hrtf_sofa, ...
    'Binauralization',          true,...
    'Downmixing_enabled',       true,...
    'Loudspeaker_setup',        struct('Shape','circular','R',2,'N',64),...
    'Rendering_mode',           'NFC_HOA',...
    'Renderer_setup',           [],...
    'Loudspeaker_type',         struct('Shape','point_source','R',0.04),...
    'Virtual_source_type',      struct('Shape','point_source','R',0.01));

if isempty(handles.sound_scene_setup.Renderer_setup)
    handles.sound_scene_setup.Renderer_setup = get_default_renderer_setup(handles.sound_scene_setup.Rendering_mode);
end
handles.sound_scene_gui = listener_space_axes(handles.axes1);
handles.sound_scene = sound_scene(handles.sound_scene_gui,handles.sound_scene_setup);

handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = Binauralizer_OutputFcn(hObject, eventdata, handles)
warning('off','all')
varargout{1} = handles.output;

% --- Executes on button press in play_btn.

function play_btn_Callback(hObject, eventdata, handles)
handles.stop_now = 0;
guidata(hObject,handles);
while (~isDone(handles.Input_stream))&&(~handles.stop_now)
    output = handles.Volume*handles.sound_scene.render_sound_scene(handles.Input_stream()...
        , handles.sound_scene_setup.Binauralization, handles.sound_scene_setup.Downmixing_enabled );
    handles.Output_device(output);
    drawnow limitrate
    handles = guidata(hObject);
end
release(handles.Output_device)
mean(diff(t))
release(handles.Input_stream)

% --- Executes on button press in load_file_btn.
function load_file_btn_Callback(hObject, eventdata, handles)
[file,path] = uigetfile('SoundSamples/*.wav;*.mp3;*.aac;*.ac3');
if file==0
    return
end
handles.Input_stream = dsp.AudioFileReader(strcat(path,file),...
    'SamplesPerFrame',handles.sound_scene_setup.Block_size);
Nch_in = handles.Input_stream.info.NumChannels;
fs = handles.Input_stream.SampleRate;
handles.sound_scene_setup.N_in = Nch_in;
handles.sound_scene_setup.SampleRate = fs;
handles.sound_scene.delete(handles.sound_scene_gui);
handles.sound_scene = sound_scene(handles.sound_scene_gui,handles.sound_scene_setup);
guidata(hObject,handles);

% --- Executes on slider movement.
function Volume_Callback(hObject, eventdata, handles)
handles.Volume = get(hObject,'Value');
guidata(hObject,handles)

function Volume_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function stop_btn_Callback(hObject, eventdata, handles)
handles.stop_now = 1;
guidata(hObject, handles);

% --- Executes on button press in harmonic_simulation.
function harmonic_simulation_Callback(hObject, eventdata, handles)
handles.simulator = [];
handles.simulator = sound_scene_simulator(handles.sound_scene, handles.sound_scene_gui, 'harmonic', 2e3 );
handles.simulator.simulate(handles.sim_t.Value);
guidata(hObject,handles);

% --- Executes on button press in impulsive_simulation.
function impulsive_simulation_Callback(hObject, ~, handles)
handles.simulator = [];
handles.simulator = sound_scene_simulator(handles.sound_scene, handles.sound_scene_gui, 'impulse', 4 );
handles.simulator.simulate(handles.sim_t.Value);
guidata(hObject,handles);

% --- Executes on slider movement.
function sim_t_Callback(hObject, eventdata, handles)
if isfield(handles,'simulator')
    handles.simulator.simulate(handles.sim_t.Value);
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sim_t_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in load_hrtf.
function load_hrtf_Callback(hObject, eventdata, handles)
% hObject    handle to load_hrtf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile(fullfile('HRTFs','*.sofa'));
if file==0
    return
end
handles.sound_scene_setup.sofa_def_path = path;
hrtf_sofa = SOFAload(fullfile(path,file));
handles.sound_scene_setup.HRTF = hrtf_sofa;
handles.sound_scene.delete(handles.sound_scene_gui);
handles.sound_scene = sound_scene(handles.sound_scene_gui,handles.sound_scene_setup);
guidata(hObject, handles);
%handles.sound_scene.update_hrtf(hrtf_sofa);
%TODO: https://www.mathworks.com/matlabcentral/answers/217751-keep-gui-functions-running-when-opening-an-uigetfile-dialog


% --- Executes on button press in Binauralization.
function Binauralization_Callback(hObject, eventdata, handles)
% hObject    handle to Binauralization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sound_scene_setup.Binauralization = logical(get(hObject,'Value'));
guidata(hObject, handles);


% --- Executes on button press in plot_transfer.
function plot_transfer_Callback(hObject, eventdata, handles)
handles.simulator = [];
handles.simulator = sound_scene_simulator(handles.sound_scene, handles.sound_scene_gui, 'impulse', 1 );
handles.simulator.plot_transfer;
guidata(hObject,handles);

