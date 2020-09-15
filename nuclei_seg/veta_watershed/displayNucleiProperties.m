function varargout = displayNucleiProperties(varargin)
% DISPLAYNUCLEIPROPERTIES MATLAB code for displayNucleiProperties.fig
%      DISPLAYNUCLEIPROPERTIES, by itself, creates a new DISPLAYNUCLEIPROPERTIES or raises the existing
%      singleton*.
%
%      H = DISPLAYNUCLEIPROPERTIES returns the handle to a new DISPLAYNUCLEIPROPERTIES or the handle to
%      the existing singleton*.
%
%      DISPLAYNUCLEIPROPERTIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISPLAYNUCLEIPROPERTIES.M with the given input arguments.
%
%      DISPLAYNUCLEIPROPERTIES('Property','Value',...) creates a new DISPLAYNUCLEIPROPERTIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before displayNucleiProperties_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to displayNucleiProperties_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help displayNucleiProperties

% Last Modified by GUIDE v2.5 09-May-2013 17:46:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @displayNucleiProperties_OpeningFcn, ...
                   'gui_OutputFcn',  @displayNucleiProperties_OutputFcn, ...
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


% --- Executes just before displayNucleiProperties is made visible.
function displayNucleiProperties_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to displayNucleiProperties (see VARARGIN)

% Choose default command line output for displayNucleiProperties
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes displayNucleiProperties wait for user response (see UIRESUME)
% uiwait(handles.mainFigure);


if length(varargin) == 4
    
    I = varargin{1};
    nuclei = varargin{2};
    properties = varargin{3};
    property = varargin{4};
    
    handles.image = image(I, 'parent', handles.imageAxes);
    axis image;
    axis off;    
    hold on;
    
    for k = 1:length(nuclei)
        
        plot(nuclei{k}(:,2), nuclei{k}(:,1), 'g-', 'LineWidth', 2);
        
    end
    
    handles.objectCentroids = cat(1,properties.Centroid);
    handles.propery = cat(1, properties.(property));
    
    guidata(handles.mainFigure, handles);
    
else
    
    error('Invalid number of inpur arguments.');
    
end

% --- Outputs from this function are returned to the command line.
function varargout = displayNucleiProperties_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on mouse motion over figure - except title and menu.
function mainFigure_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to mainFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentPointInFigure = get(hObject, 'CurrentPoint');

axesPosition = get(handles.imageAxes, 'Position');

isInAxes = (currentPointInFigure(1,1) > axesPosition(1)) && ...
    (currentPointInFigure(1,1) < axesPosition(1)+axesPosition(3)) && ...
    (currentPointInFigure(1,2) > axesPosition(2)) && ...
    (currentPointInFigure(1,2) < axesPosition(2)+axesPosition(4));

if isInAxes
    
    currentPointInAxes = get(handles.imageAxes, 'CurrentPoint');
    
    XData = get(handles.image, 'XData');
    YData = get(handles.image, 'YData');
    
    isInImage = (currentPointInAxes(1,1) > 0) && ...
        (currentPointInAxes(1,1) < XData(2)) && ...
        (currentPointInAxes(1,2) > 0) && ...
        (currentPointInAxes(1,2) < YData(2));
    
    if isInImage
        
        idx = findNearestObject(handles.objectCentroids, currentPointInAxes(1,1:2));
        
        currentString = num2str(handles.propery(idx), '%2.2f');

        set(handles.propertyDisplayEdit, 'String', currentString);
        
    end
    
end

% --- Executes during object creation, after setting all properties.
function propertyDisplayEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to propertyDisplayEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function idx = findNearestObject(objectCentroids, currentPoint)

D = bsxfun(@minus, objectCentroids, currentPoint);
D = sum(abs(D), 2);

[~, idx] = min(D);
   
