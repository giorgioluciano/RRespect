function varargout = specgui(varargin)
% SPECGUI MATLAB code for specgui.fig
%      SPECGUI, by itself, creates a new SPECGUI or raises the existing
%      singleton*.
%
%      H = SPECGUI returns the handle to a new SPECGUI or the handle to
%      the existing singleton*.
%
%      SPECGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPECGUI.M with the given input arguments.
%
%      SPECGUI('Property','Value',...) creates a new SPECGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before specgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to specgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help specgui

% Last Modified by GUIDE v2.5 20-Apr-2012 12:44:05
addpath('./gui:./common:./continuous:./disc:./output')
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @specgui_OpeningFcn, ...
                   'gui_OutputFcn',  @specgui_OutputFcn, ...
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


% --- Executes just before specgui is made visible.
function specgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to specgui (see VARARGIN)

handles.CalcContSpec = 0;
handles.CalcDiscSpec = 0;

% Choose default command line output for specgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes specgui wait for user response (see UIRESUME)
% uiwait(handles.guimain);


% --- Outputs from this function are returned to the command line.
function varargout = specgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in buttoncompute.
function buttoncompute_Callback(hObject, eventdata, handles)
% hObject    handle to buttoncompute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%
% Continuous or discrete spectrum
%

IsContSpec=get(handles.radiobuttonCont,'Value');

%
% Standard settings regardless of choice
%
par.verbose  = 1;
par.plotting = 0;
par.GstFile  = 'Gst.dat';

%
% Continuous
%
if IsContSpec
    
    par.ns = str2num(get(handles.editns,'String'));
    isauto = get(handles.buttoncontauto,'Value');

    %
    % Choose lambdaC automatically
    %
    if isauto
        
        par.lamC = 0;
        par.SmFacLam = str2num(get(handles.textSF,'String'));
        
    else  % read supplied lambda
        
        par.lamC = str2num(get(handles.editlamC,'String'));
        if isnan(par.lamC) || par.lamC < 0
            errordlg('lambdaC needs to be specified. Running in auto mode');
            par.lamC = 0; 
            par.SmFacLam = str2num(get(handles.textSF,'String'));
        end
        par.SmFacLam = 0;
        
    end

    par.FreqEnd = 1;
    
    %
    % Set up actual calculation. Update trigger
    %
    set(handles.textstatus,'String','Computing!');
    contSpec(par);
    set(handles.textstatus,'String','Continuous Spectrum Computed');
    handles.CalcContSpec = 1;

    updateplot(handles.popupmenugraph1,handles.axes1);
    updateplot(handles.popupmenugraph2,handles.axes2);

else  % discrete spectrum
    
    if(handles.CalcContSpec) % only if continuous spectrum computed before
        
        isauto = get(handles.buttondiscauto,'Value');

        if isauto
            par.Nopt = 0;
        else
            par.Nopt = str2num(get(handles.editNopt,'String'));
        end
        
        if(get(handles.buttonprune,'Value'))
            par.prune = 1;
        else
            par.prune = 0;
        end
        
        par.BaseDistWt = str2num(get(handles.textflat,'String'));
        par.condWt     = str2num(get(handles.textcond,'String'));
        discSpec(par);
        set(handles.textstatus,'String','Continuous and Discrete Spectrum Computed');
        handles.CalcDiscSpec = 1; 
        
    else
        errordlg('Need to Compute Continuous Spectrum First')
    end
    
    updateplot(handles.popupmenugraph1,handles.axes1);
    updateplot(handles.popupmenugraph2,handles.axes2);
    
end


guidata(hObject,handles);

% --- Executes on slider movement.
function sliderSF_Callback(hObject, eventdata, handles)
% hObject    handle to sliderSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
num=get(handles.sliderSF,'Value');
set(handles.textSF,'String',num2str(num))
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderSF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editns_Callback(hObject, eventdata, handles)
% hObject    handle to editns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editns as text
%        str2double(get(hObject,'String')) returns contents of editns as a double


% --- Executes during object creation, after setting all properties.
function editns_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in buttoncontauto.
function buttoncontauto_Callback(hObject, eventdata, handles)
% hObject    handle to buttoncontauto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buttoncontauto


% --- Executes on button press in buttoncontmanual.
function buttoncontmanual_Callback(hObject, eventdata, handles)
% hObject    handle to buttoncontmanual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buttoncontmanual


% --- Executes on slider movement.
function sliderflat_Callback(hObject, eventdata, handles)
% hObject    handle to sliderflat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
num=get(handles.sliderflat,'Value');
set(handles.textflat,'String',num2str(num))
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderflat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderflat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editNopt_Callback(hObject, eventdata, handles)
% hObject    handle to editNopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNopt as text
%        str2double(get(hObject,'String')) returns contents of editNopt as a double


% --- Executes during object creation, after setting all properties.
function editNopt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slidercond_Callback(hObject, eventdata, handles)
% hObject    handle to slidercond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
num=get(handles.slidercond,'Value');
set(handles.textcond,'String',num2str(num))
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slidercond_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slidercond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editlamC_Callback(hObject, eventdata, handles)
% hObject    handle to editlamC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editlamC as text
%        str2double(get(hObject,'String')) returns contents of editlamC as a double


% --- Executes during object creation, after setting all properties.
function editlamC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editlamC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenugraph1.
function popupmenugraph1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenugraph1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenugraph1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenugraph1
ichoice = get(hObject,'Value');
if ichoice <= 3 && handles.CalcContSpec
    updateplot(hObject,handles.axes1)
elseif ichoice > 3 && handles.CalcDiscSpec
    updateplot(hObject,handles.axes1)
end

% --- Executes during object creation, after setting all properties.
function popupmenugraph1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenugraph1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenugraph2.
function popupmenugraph2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenugraph2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenugraph2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenugraph2
ichoice = get(hObject,'Value');
if ichoice <= 3 && handles.CalcContSpec
    updateplot(hObject,handles.axes2)
elseif ichoice > 3 && handles.CalcDiscSpec
    updateplot(hObject,handles.axes2)
end

% --- Executes during object creation, after setting all properties.
function popupmenugraph2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenugraph2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function textSF_Callback(hObject, eventdata, handles)
% hObject    handle to textSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textSF as text
%        str2double(get(hObject,'String')) returns contents of textSF as a double
num=str2num(get(handles.textSF,'String'));
if (isnan(num) || num > 1 || num < -1)
    num = 0.0;
    set(hObject,'String',num);
    errordlg('Input must be between -1 and 1', 'Error')
else
    set(handles.sliderSF,'Value',num);    
end
guidata(hObject,handles) 

% --- Executes during object creation, after setting all properties.
function textSF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function textflat_Callback(hObject, eventdata, handles)
% hObject    handle to textflat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textflat as text
%        str2double(get(hObject,'String')) returns contents of textflat as a double
num=str2num(get(handles.textflat,'String'));
if (isnan(num) || num > 1 || num < 0)
    num = 0.5;
    set(hObject,'String',num);
    errordlg('Input must be between 0 and 1', 'Error')
else
    set(handles.sliderflat,'Value',num);
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function textflat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textflat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function textcond_Callback(hObject, eventdata, handles)
% hObject    handle to textcond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of textcond as text
%        str2double(get(hObject,'String')) returns contents of textcond as a double
num=str2num(get(handles.textcond,'String'));
if (isnan(num) || num > 1 || num < 0)
    num = 0.5;
    set(hObject,'String',num);
    errordlg('Input must be between 0 and 1', 'Error')
else
    set(handles.slidercond,'Value',num)
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function textcond_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textcond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
