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
addpath(genpath('HRTFs'))
SOFAstart;
hrtf_sofa = SOFAload(['BuK_ED_corr.sofa']);

% Setup options:
%   Rendering:              'Binaural':        Binaural rendering only
%                           'WFS':             Wave Field Synthesis driving
%                                              functions
%                           'VBAP':            Vector-based Amplitude
%                                              Panning
%                           '':
%   Binaural_source_type:   'point_source':    omnidirectional directivity,
%                                              parameter: R (only for drawing)
%                           'circular_piston': baffled piston directivity,
%                                              parameter: R radius
%                           'two_way_speaker': two baffled pistons, with
%                                              Linkwitz-Riley crossover filter,
%                                              parameter: [R_lp, R_hp] radii
%   Virtual_source_type:    'point_source':    omnidirectional directivity,
%                                              parameter: R (only for drawing)
%                           'plane_wave':      plane_,
%                                              parameter: R radius

N_ssd = 48;
dx = 2*2*pi/N_ssd;
r0  = 2.22/pi*dx;
if isempty(varargin)
    input_file = 'gitL.wav';
else
    input_file = varargin{1};
end
handles.sound_scene_setup = struct(  ...
    'Input_file',               input_file,...
    'Block_size',               1024, ...
    'HRTF',                     hrtf_sofa, ...
    'sofa_def_path',            '..',...
    'Volume',                   0.5, ...
    'Binaural_source_type',     struct('Shape','circular_piston','R',0.06),...
    'Virtual_source_type',      struct('Shape','point_source','R',0.01),...
    'Rendering',                'WFS',...
    'renderer_setup',           struct('R',2,'N',N_ssd,'Antialiasing','off'));

if length(varargin) ~= 0
    handles.sound_scene_setup.Rendering = 'Binaural';
end
handles.sound_scene_setup.Input_stream = dsp.AudioFileReader(handles.sound_scene_setup.Input_file,...
    'SamplesPerFrame',handles.sound_scene_setup.Block_size,'PlayCount',10);
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
deviceWriter = audioDeviceWriter('SampleRate',handles.sound_scene_setup.Input_stream.SampleRate);
guidata(hObject,handles);
elapsed_time = 0;
i = 1;
while (~isDone(handles.sound_scene_setup.Input_stream))&&(~handles.stop_now)
    tic
    if handles.Bypass.Value == 1
        output = handles.sound_scene_setup.Input_stream();
    else
    output = handles.sound_scene.binauralize_sound_scene(handles.sound_scene_setup.Volume*...
        handles.sound_scene_setup.Input_stream());
    end
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
if file==0
  return
end
handles.sound_scene_setup.Input_file = strcat(path,file);
handles.sound_scene_setup.Input_stream = dsp.AudioFileReader(handles.sound_scene_setup.Input_file,...
    'SamplesPerFrame',handles.sound_scene_setup.Block_size);
handles.sound_scene.delete(handles.sound_scene_gui);
handles.sound_scene = sound_scene(handles.sound_scene_gui,handles.sound_scene_setup);
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

% --- Executes on button press in harmonic_simulation.
function harmonic_simulation_Callback(hObject, eventdata, handles)
handles.simulator = [];
handles.simulator = sound_scene_simulator(handles.sound_scene, handles.sound_scene_gui, 'harmonic', 1e3 );
handles.simulator.simulate(handles.sim_t.Value);
guidata(hObject,handles);

% --- Executes on button press in impulsive_simulation.
function impulsive_simulation_Callback(hObject, eventdata, handles)
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
[file,path] = uigetfile(fullfile(handles.sound_scene_setup.sofa_def_path,'*.sofa'));
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


% --- Executes on button press in Bypass.
function Bypass_Callback(hObject, eventdata, handles)