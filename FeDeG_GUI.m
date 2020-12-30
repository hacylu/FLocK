function varargout = FeDeG_GUI(varargin)
% FEDEG_GUI MATLAB code for FeDeG_GUI.fig
%      FEDEG_GUI, by itself, creates a new FEDEG_GUI or raises the existing
%      singleton*.
%
%      H = FEDEG_GUI returns the handle to a new FEDEG_GUI or the handle to
%      the existing singleton*.
%
%      FEDEG_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FEDEG_GUI.M with the given input arguments.
%
%      FEDEG_GUI('Property','Value',...) creates a new FEDEG_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FeDeG_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FeDeG_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FeDeG_GUI

% Last Modified by GUIDE v2.5 19-Jul-2020 10:34:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FeDeG_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @FeDeG_GUI_OutputFcn, ...
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


% --- Executes just before FeDeG_GUI is made visible.
function FeDeG_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FeDeG_GUI (see VARARGIN)

% Choose default command line output for FeDeG_GUI
handles.output = hObject;
if ispc()
    addpath(genpath([pwd '\nuclei_seg']));
else
    addpath(genpath([pwd '/nuclei_seg']));
end
handles.imgfilename = [];
handles.imgdata = [];
handles.nuclei = [];
handles.properties = [];
handles.feature_space= 'Centroid-Area';
handles.bandWidth_space=20;
handles.bandWidth_features=20;
handles.num_fixed_types = 0;

handles.debug=[]; 
handles.shownuclei_thiness=[]; 
handles.shownucleiconnection_thiness=[]; 

handles.data_other_attribute = [];
handles.clust2types = [];
handles.typeCent = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FeDeG_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FeDeG_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% read image
msgbox('Please select appropriate size within 1MB of the image,otherwise it will take a long time to process the image');
[imgfilename,pathname]=uigetfile({'*.bmp;*.jpg;*.png;*.jpeg;*.tif'},'select a image');
str=[pathname imgfilename];
 if imgfilename
    imgdata = imread(str);
    figure('name','Original image'); 
    imshow(imgdata); 
    handles.imgfilename = imgfilename;
    handles.imgdata = imgdata;
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pushbutton1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% nuclei segmentation
if ~isempty(handles.imgdata)
    [I_norm,~,~] = normalizeStaining(handles.imgdata);
    I_normRed=I_norm(:,:,1);
    p.scales=3:2:10;
    [nuclei, properties] = nucleiSegmentationV2(single(I_normRed),p);  % 灰度图像要进行类型转换
    figure('name','Segmentation result');
    imshow(handles.imgdata);
    hold on;
    for k = 1:length(nuclei)
        plot(nuclei{k}(:,2), nuclei{k}(:,1), 'g-', 'LineWidth', 2);
    end
    hold off;
    handles.nuclei = nuclei;
    handles.properties = properties;
end
guidata(hObject, handles);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% construct FeDeG and extract feature
handles.debug=1; 
handles.shownuclei_thiness=1; 
handles.shownucleiconnection_thiness=3; 

if ~isempty(handles.feature_space) && ~isempty(handles.bandWidth_space) && ~isempty(handles.bandWidth_features) &&~isempty(handles.num_fixed_types)
    disp('use MS to build FeDeG...');
    [clustCent,data2cluster,cluster2dataCell,data_other_attribute,clust2types,typeCent] = Lconstruct_FeDeG_v2(handles.imgdata,handles);

    handles.imgdata = handles.imgdata;
    handles.debug = 1;
    handles.data_other_attribute = data_other_attribute;
    handles.clust2types = clust2types;
    handles.typeCent = typeCent;
    [set_feature,set_feature_name] = L_get_FeDeG_features_v2(clustCent,cluster2dataCell,handles);
    
    % computer typical set_feature_name's value
    idx1 = strcmp(set_feature_name,'FeDeG-portion of intersected FeDeG');
    str1 = ['FeDeG-portion of intersected FeDeG:  ',num2str(set_feature(idx1))];

    idx2 = strcmp(set_feature_name,'FeDeG-abs. number of intersected FeDeG');
    str2 = ['FeDeG-abs. number of intersected FeDeG:  ',num2str(set_feature(idx2))];
    
    idx3 = strcmp(set_feature_name,'FeDeG-portion of highly intersected FeDeG:T=0.3');
    str3 = ['FeDeG-portion of highly intersected FeDeG:T=0.3:  ',num2str(set_feature(idx3))];
    
    idx4 = strcmp(set_feature_name,'FeDeG-mean(size of non-1-2cell FeDeG)');
    str4 = ['FeDeG-mean(size of non-1-2cell FeDeG):  ',num2str(set_feature(idx4))];
    
    idx5 = strcmp(set_feature_name,'FeDeG-min(nuclei no. in FeDeG)');
    str5 = ['FeDeG-min(nuclei no. in FeDeG):  ',num2str(set_feature(idx5))];
    
    idx6 = strcmp(set_feature_name,'FeDeG-max(nuclei no. in FeDeG)');
    str6 = ['FeDeG-max(nuclei no. in FeDeG):  ',num2str(set_feature(idx6))];
    str7 = []; str8 = [];
    if handles.num_fixed_types == 2
        idx7 = strcmp(set_feature_name,'FeDeG-mean(var of nuclear feat1to the centroid)');
        str7 = ['FeDeG-mean(var of nuclear feat1to the centroid):  ',num2str(set_feature(idx7))];

        idx8 = strcmp(set_feature_name,'FeDeG-mean(var of nuclear feat2to the centroid)');
        str8 = ['FeDeG-mean(var of nuclear feat2to the centroid):  ',num2str(set_feature(idx8))];
    end
    %str = [str1 str2 str3 str4 str5 str6 str7 str8];
    msgbox({'Selected feature result:',' ',str1,str2,str3,str4,str5,str6,str7,str8},'Selected feature result','help');
end
guidata(hObject, handles);

% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get feature_space
str = get(handles.popupmenu3,'String');
var = get(handles.popupmenu3,'value');

switch var
    case 1
        handles.feature_space = str{1};
    case 2
        handles.feature_space=str{2};
    case 3
        handles.feature_space=str{3};
    case 4
        handles.feature_space=str{4};
    case 5
        handles.feature_space=str{5};
    case 6
        handles.feature_space=str{6};
    case 7
        handles.feature_space=str{7};
    case 8
        handles.feature_space=str{8};
end
guidata(hObject, handles);


% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get bandWidth_space
str1 = get(handles.popupmenu4,'String');
var1 = get(handles.popupmenu4,'value');
switch var1
    case 1
        handles.bandWidth_space=str2num(str1{1});
    case 2
        handles.bandWidth_space=str2num(str1{2});
    case 3
        handles.bandWidth_space=str2num(str1{3});
    case 4
        handles.bandWidth_space=str2num(str1{4});
    case 5
        handles.bandWidth_space=str2num(str1{5});
end
guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get bandWidth_features
var = get(handles.popupmenu3,'value');

switch var
    case 1
        answer = inputdlg('area ([10,100])','Please input...');
        handles.bandWidth_features = str2num(answer{1});
        set(handles.text26,'string',['area = ',num2str(handles.bandWidth_features)]);
    case 2
        answer = inputdlg({'area ([10,100])','Eccentricity([0.05,0.2])'},'Please input...');
        handles.bandWidth_features = [str2num(answer{1});str2num(answer{2})];
        set(handles.text26,'string',['area = ',num2str(handles.bandWidth_features(1)),' ; ','Eccentricity = ',num2str(handles.bandWidth_features(2))]);
    case 3
        answer = inputdlg({'area ([10,100])','MeanIntensity ([5,80])'},'Please input...');
        handles.bandWidth_features = [str2num(answer{1});str2num(answer{2})];
        set(handles.text26,'string',['area = ',num2str(handles.bandWidth_features(1)),' ; ','MeanIntensity = ',num2str(handles.bandWidth_features(2))]);
    case 4
        answer = inputdlg('MeanIntensity ([5,80])','Please input...');
        handles.bandWidth_features = str2num(answer{1});
        set(handles.text26,'string',['MeanIntensity = ',num2str(handles.bandWidth_features)]);
    case 5
        answer = inputdlg('Eccentricity ([0.05,0.2])','Please input...');
        handles.bandWidth_features = str2num(answer{1});
        set(handles.text26,'string',['Eccentricity = ',num2str(handles.bandWidth_features)]);
    case 6
        answer = inputdlg('Longness ([0.25,2])','Please input...');
        handles.bandWidth_features = str2num(answer{1});
        set(handles.text26,'string',['Longness = ',num2str(handles.bandWidth_features)]);
    case 7
        answer = inputdlg('Circularity ([0.05,0.2])','Please input...');
        handles.bandWidth_features = str2num(answer{1});
        set(handles.text26,'string',['Circularity = ',num2str(handles.bandWidth_features)]);
    case 8
        answer = inputdlg({'area ([10,100])','Eccentricity ([0.05,0.2])','MeanIntensity ([5,80])'},'Please input...');
        handles.bandWidth_features = [str2num(answer{1});str2num(answer{2});str2num(answer{3})];
        set(handles.text26,'string',['area = ',num2str(handles.bandWidth_features(1)),' ; ','Eccentricity = ',num2str(handles.bandWidth_features(2)),' ; ','MeanIntensity = ',num2str(handles.bandWidth_features(3))]);
end
guidata(hObject, handles);


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get num_fixed_types
answer = inputdlg('specify the phenotype numbers in the image','please input...');
handles.num_fixed_types = str2num(answer{1});
set(handles.text27,'string',['num_fixed_types = ',num2str(handles.num_fixed_types)]);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function text26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
