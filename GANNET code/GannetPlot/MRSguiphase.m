function varargout = MRSguiphase(varargin)
% MRSGUIPHASE M-file for MRSguiphase.fig
%      MRSGUIPHASE, by itself, creates a new MRSGUIPHASE or raises the existing
%      singleton*.
%
%      H = MRSGUIPHASE returns the handle to a new MRSGUIPHASE or the handle to
%      the existing singleton*.
%
%      MRSGUIPHASE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MRSGUIPHASE.M with the given input arguments.
%
%      MRSGUIPHASE('Property','Value',...) creates a new MRSGUIPHASE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MRSguiphase_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MRSguiphase_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% CJE 110215: First version.  Now stable

% Edit the above text to modify the response to help MRSguiphase

% Last Modified by GUIDE v2.5 13-Apr-2011 22:33:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRSguiphase_OpeningFcn, ...
                   'gui_OutputFcn',  @MRSguiphase_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before MRSguiphase is made visible.
function MRSguiphase_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MRSguiphase (see VARARGIN)

% Choose default command line output for MRSguiphase
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(hObject, 'Name', 'GannetPlot');

% Populate listboxes and phases
update_listbox1(handles);
update_listbox2(handles);
set_phase_from_struct(handles);
plot_all_spectra(handles);

% UIWAIT makes MRSguiphase wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MRSguiphase_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;





%%%%%%%%%%%%%%%%%%%% WORKSPACE STRUCT LISTBOX %%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
update_listbox1(handles)
set(handles.listbox2,'Value',1); % set to first p file in struct
update_listbox2(handles)
set_phase_from_struct(handles)
plot_one_spectrum(handles)


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');

%update_listbox;

end


%%% CJE %%%
function update_listbox1(handles)
allvars = evalin('base','who');
numstructs=0;
for jj = 1 : length(allvars)
  % char(39) is single quote '
  cmd1 = ['isfield(' allvars{jj} ', ' char(39) 'pfile' char(39) ')'];
  structisMRS = evalin('base', cmd1);
  if(structisMRS)
    numstructs=numstructs+1;
    MRSstructlist{numstructs}=allvars{jj};
  end
end
if numstructs == 0
  MRSstructlist{1} = 'No_Data';
end

set(handles.listbox1,'String',MRSstructlist);

function [var1,var2] = get_pfile_list(handles)
% Returns the names of the two variables to plot
list_entries = get(handles.listbox1,'String');
%index_selected = get(handles.listbox1,'Value');



%%%%%%%%%%%%%%%%% PFILE STRUCT LISTBOX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2
update_listbox1(handles)
update_listbox2(handles)
set_phase_from_struct(handles)
plot_one_spectrum(handles)

% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%  update listbox2 - find pfiles in struct %%%%
function update_listbox2(handles)
list_entries = get(handles.listbox1,'String');
index_selected = get(handles.listbox1,'Value');


nameofstruct=list_entries{index_selected};
nameofstructpfiles=[nameofstruct '.pfile'];

% catch 'No Data' condition
if (strcmp(nameofstruct,'No_Data'))
  structisMRS = 0;
else
  cmd1 = ['isfield(' nameofstruct ', ' char(39) 'pfile' char(39) ')'];
  structisMRS = evalin('base', cmd1);
end

if(structisMRS) 
    nameofstructpfiles=[nameofstruct '.pfile'];
    text11=evalin('base', nameofstructpfiles);
    set(handles.listbox2,'String',text11)
else set(handles.listbox2,'String','Not MRS data')
end


%%%%%%%%%%%%%%%%%%%%%%% PLOT ALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_all_spectra(handles)



%%%%%%%%%%%%%%%%%%%%% PLOT ONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Plot One
plot_one_spectrum(handles)

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_sum_spectra(handles)


%%%%%%%%%%%%%%%%%%%%%%% ZEROTH PHASE SLIDER %%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on slider movement.
function slider_phase0_Callback(hObject, eventdata, handles)
% hObject    handle to slider_phase0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
update_edit_phase0(handles);
phase_spectrum(handles);


% --- Executes during object creation, after setting all properties.
function slider_phase0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_phase0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% update the slider phase based on value typed into textbox
function update_slider_phase0(handles)
phase0str=get(handles.edit_phase0, 'String');
phase0val=str2num(phase0str);
if isempty(phase0val)
  % not a float/int - don't update
  update_edit_phase0(handles)
else
  set(handles.slider_phase0, 'Value',phase0val);
  phase0str2 = sprintf('%.1f', phase0val);
  set(handles.edit_phase0, 'String', phase0str2);
end


%%%%%%%%%%%%%%%%%%%% FIRST ORDER PHASE SLIDER %%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on slider movement.
function slider_phase1_Callback(hObject, eventdata, handles)
% hObject    handle to slider_phase1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
update_edit_phase1(handles);
phase_spectrum(handles);


% --- Executes during object creation, after setting all properties.
function slider_phase1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_phase1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function update_slider_phase1(handles)
phase1str=get(handles.edit_phase1, 'String');
phase1val=str2num(phase1str);
if isempty(phase1val)
  % not a float/int - don't update
  update_edit_phase1(handles)
else
  set(handles.slider_phase1, 'Value',phase1val);
  phase1str2 = sprintf('%.1f', phase1val);
  set(handles.edit_phase1, 'String', phase1str2);
end


%%%%%%%%%%%%%%%%%%%% ZEROTH PHASE EDIT BOX %%%%%%%%%%%%%%%%%%%%
 % --- Executes during object creation, after setting all properties.
function edit_phase0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_phase0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
  
function edit_phase0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_phase0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_phase0 as text
%        str2double(get(hObject,'String')) returns contents of
%        edit_phase0 as a double

%Also need to update the slider based on values typed into phase
%box

update_slider_phase0(handles);
phase_spectrum(handles);

%update textbox phase val based on slider
function update_edit_phase0(handles)
phase0val=get(handles.slider_phase0, 'Value');
phase0str=sprintf('%.1f', phase0val);
set(handles.edit_phase0, 'String', phase0str);


%%%%%%%%%%%%%%%%%%%% FIRST ORDER PHASE EDIT BOX %%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edit_phase1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_phase1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_phase1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_phase1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_phase1 as text
%        str2double(get(hObject,'String')) returns contents of edit_phase1 as a double

%Also need to update the slider based on values typed into phase
%box

update_slider_phase1(handles);
phase_spectrum(handles);


% --- Executes on button press in LoadPfile.
function LoadPfile_Callback(hObject, eventdata, handles)
% hObject    handle to LoadPfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ffilter = { '*.7' 'GE P-file'};
[fname, fpath] = uigetfile(ffilter, 'Open P-file');
if isequal(fname,0) | isequal(fpath,0)
  %disp('')
else
  fullpath = [fpath fname];
  pos = regexp(fname, '\.7')
  structname = ['pfile_' fname((pos-6):(pos-1))];
  %structname = ['pfile'];
  cmd = [structname '= MRSLoadPfiles({' char(39) fullpath char(39) '},0);'  ]
  evalin('base',cmd);
  
  update_listbox1(handles);
  update_listbox2(handles);
  set_phase_from_struct(handles);
  plot_one_spectrum(handles);

end


%%%%%%%%%%%%%%%%%%%%%%%%% CJE functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_edit_phase1(handles)
phase1val=get(handles.slider_phase1, 'Value');
phase1str=sprintf('%.1f', phase1val);
set(handles.edit_phase1, 'String', phase1str);

function plot_all_spectra(handles)
% Plot All
list_entries = get(handles.listbox1,'String');
index_selected = get(handles.listbox1,'Value');

if (strcmp(list_entries{index_selected},'No_Data')==0)
  nameofstruct=list_entries{index_selected};
  nameofstructpfiles=[nameofstruct '.pfile'];
  cmd1 = ['isfield(' nameofstruct ', ' char(39) 'pfile' char(39) ')'];
  structisMRS = evalin('base', cmd1);

  if(structisMRS) 
    axes(handles.axes_MRS)
    title('hello')
    cmd2 = ['MRSplotstack(' list_entries{index_selected} ')' ];
    evalin('base',cmd2);
  end
end

function plot_one_spectrum(handles)
% Plot All
list_entries = get(handles.listbox1,'String');
index_selected = get(handles.listbox1,'Value');

if (strcmp(list_entries{index_selected},'No_Data')==0)
  spec_selected= get(handles.listbox2,'Value');
  nameofstruct=list_entries{index_selected};
  nameofstructpfiles=[nameofstruct '.pfile'];
  cmd1 = ['isfield(' nameofstruct ', ' char(39) 'pfile' char(39) ')'];
  structisMRS = evalin('base', cmd1);

  if(structisMRS) 
    axes(handles.axes_MRS)
    cmd2 = ['MRSplotstackone(' list_entries{index_selected} ',' ...
	    num2str(spec_selected) ')' ] ;
    evalin('base',cmd2);
  end
end


function plot_sum_spectra(handles)
% Plot All
list_entries = get(handles.listbox1,'String');
index_selected = get(handles.listbox1,'Value');

if (strcmp(list_entries{index_selected},'No_Data')==0)  
  nameofstruct=list_entries{index_selected};
  nameofstructpfiles=[nameofstruct '.pfile'];
  cmd1 = ['isfield(' nameofstruct ', ' char(39) 'pfile' char(39) ')'];
  structisMRS = evalin('base', cmd1);

  if(structisMRS) 
    axes(handles.axes_MRS)
    title('hello')
    cmd2 = ['MRSplotstacksum(' list_entries{index_selected} ')' ];
    evalin('base',cmd2);
  end

end


% phase the spectrum, based on values in sliders and currently
% highlighted P file
function phase_spectrum(handles)
%index of highlighted P file
list_entries = get(handles.listbox1,'String');
index_selected = get(handles.listbox1,'Value');

if (strcmp(list_entries{index_selected},'No_Data')==0)
  activestruct =list_entries{index_selected};
  activespec = get(handles.listbox2,'Value'); 
  phase0val=get(handles.slider_phase0, 'Value');
  phase1val=get(handles.slider_phase1, 'Value');

  cmd2 = [activestruct '= MRSphase_set(' ...
	  activestruct ',' ... 
	  num2str(activespec) ',' ...
	  num2str(phase0val) ',' ...
	  num2str(phase1val) ')' ];
  evalin('base',cmd2);

  plot_one_spectrum(handles);
end


function set_phase_from_struct(handles)
%read then refresh phase textboxes and sliders from 
% saved phase values in struct
list_entries = get(handles.listbox1,'String');
index_selected = get(handles.listbox1,'Value');
activestruct =list_entries{index_selected};
activespec = get(handles.listbox2,'Value') ;

if (strcmp(list_entries{index_selected},'No_Data')==0)
  cmd = [ activestruct '.phase(' num2str(activespec) ')' ];
  start_phase = evalin('base', cmd);

  cmd = [ activestruct '.phase_firstorder(' num2str(activespec) ')' ];
  start_phase1 = evalin('base', cmd);

  %set the textbox
  phase0str=sprintf('%.1f', start_phase);
  set(handles.edit_phase0, 'String', start_phase);
  %then update slider
  update_slider_phase0(handles);

  phase1str=sprintf('%.1f', start_phase1);
  set(handles.edit_phase1, 'String', start_phase1);
  %then update slider
  update_slider_phase1(handles);
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(handles.figure1);
% this exits matlab - use on scanner, maybe?
%exit

