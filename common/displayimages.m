function varargout = displayimages(varargin)
% DISPLAYIMAGES M-file for displayimages.fig
%      DISPLAYIMAGES, by itself, creates a new DISPLAYIMAGES or raises the existing
%      singleton*.
%
%      H = DISPLAYIMAGES returns the handle to a new DISPLAYIMAGES or the handle to
%      the existing singleton*.
%
%      DISPLAYIMAGES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISPLAYIMAGES.M with the given input arguments.
%
%      DISPLAYIMAGES('Property','Value',...) creates a new DISPLAYIMAGES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before displayimages_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to displayimages_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help displayimages

% Last Modified by GUIDE v2.5 01-Jun-2016 15:43:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @displayimages_OpeningFcn, ...
                   'gui_OutputFcn',  @displayimages_OutputFcn, ...
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


% --- Executes just before displayimages is made visible.
function displayimages_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to displayimages (see VARARGIN)

% Choose default command line output for displayimages
handles.output = hObject;
handles.figure = hObject;

handles.images1 = cell2mat(varargin(1));
handles.images2 = cell2mat(varargin(2));
titlestr = char(varargin(3));
set(handles.figure, 'Name', titlestr, 'MenuBar', 'figure');

[H,W,D] = size(handles.images1);

handles.nimages = min(size(handles.images1,3),size(handles.images2,3));
handles.imagemax = 1.25 * max( max(handles.images1(:)),max(handles.images2(:)) );
min1 = min( min(handles.images1(:)),min(handles.images2(:)) );
if (min1 < 0 )
    handles.imagemin = -1.25 * min1;
else
    handles.imagemin = 0.5 * min1;
end

sl1val = 1.0;
if (handles.nimages > 1)
    set(handles.slider1,'Min',0,'Max',1.0,...
          'Value',sl1val,'SliderStep',[(1/(handles.nimages-1)) (2/(handles.nimages-1))]);    % Image scroller
else
    set(handles.slider1,'Min',0,'Max',1.0,...
          'Value',sl1val,'SliderStep',[0.25 0.5]);    % Image scroller
end
handles.imageindex = handles.nimages;
set(handles.edit1,'Value',handles.imageindex);    % Image number update

handles.displaymax = handles.imagemax;
set(handles.edit3,'String',sprintf('%.1f',handles.displaymax) );    % displaymax update
handles.displaymin = handles.imagemin;
set(handles.edit2,'String',sprintf('%.1f',handles.displaymin) );    % displaymin update

handles.dispimage1 =  handles.images1(:,:,handles.imageindex);
handles.dispimage2 =  handles.images2(:,:,handles.imageindex);
axes(handles.axes1);
imagesc( handles.dispimage1,[handles.displaymin,handles.displaymax] ), colormap(gray), colorbar, axis square, axis off;
axes(handles.axes2);
imagesc( handles.dispimage2,[handles.displaymin,handles.displaymax] ), colormap(gray), colorbar, axis square, axis off;

img = handles.dispimage1;
maximg = max(img(:));
mask = img*0.0;
mask( img > 0.1*maximg )= 1.0;
handles.disproi = mask;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes displayimages wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = displayimages_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

sl1min = get(hObject,'Min');
sl1max = get(hObject,'Max');
sl1pos = get(hObject,'Value');
sl1val = (sl1pos/(sl1max-sl1min));

imageindex = fix(sl1val*(handles.nimages-1))+1;
if (imageindex < 1)
    imageindex = 1;
elseif (imageindex > handles.nimages)
    imageindex = handles.nimages;
end
imageindexstr = sprintf('%i',imageindex);
handles.imageindex = imageindex;
set(handles.edit1,'String',imageindexstr);    % Image number update

axes(handles.axes1);
imagesc( handles.images1(:,:,handles.imageindex),[handles.displaymin,handles.displaymax] ), colormap(gray), colorbar, axis square, axis off;
axes(handles.axes2);
imagesc( handles.images2(:,:,handles.imageindex),[handles.displaymin,handles.displaymax] ), colormap(gray), colorbar, axis square, axis off;


% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figpath = uigetdir('Choose directory to Save Figure');
filename = [figpath filesep handles.titlestr '.jpg' ];
print( '-noui', '-djpeg', filename );



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
handles.displaymin = str2double(get(hObject,'String'));
axes(handles.axes1);
imagesc( handles.images1(:,:,handles.imageindex),[handles.displaymin,handles.displaymax] ), colormap(gray), colorbar, axis square, axis off;
axes(handles.axes2);
imagesc( handles.images2(:,:,handles.imageindex),[handles.displaymin,handles.displaymax] ), colormap(gray), colorbar, axis square, axis off;

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
handles.displaymax = str2double(get(hObject,'String'));
axes(handles.axes1);
imagesc( handles.images1(:,:,handles.imageindex),[handles.displaymin,handles.displaymax] ), colormap(gray), colorbar, axis square, axis off;
axes(handles.axes2);
imagesc( handles.images2(:,:,handles.imageindex),[handles.displaymin,handles.displaymax] ), colormap(gray), colorbar, axis square, axis off;

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ( ~isfield(handles, 'disproi') )
    disp(' No dispROI found. Click dispROI first');
else
    resultsdir = uigetdir('.','Choose a directory to save excel outputfile');
    filename = input('Input excel file name :','s');
    fullfilename = fullfile(resultsdir,[filename '-regionalCESTvalues.csv']);
    fid = fopen(fullfilename,'w');
    fprintf(fid,'%s\n\n',fullfilename);
    
    fprintf(fid,'%s\n','imageindex,meanimage1,meanimage2');
    fprintf('%s\n','imageindex meanimage1 meanimage2');
    mask = handles.disproi;
    maskindex = find( mask > 0.01);
    count = length(maskindex);
    for i = 1:handles.nimages
        sum1 = 0.0;
        sum2 = 0.0;
        im1 = handles.images1(:,:,i);
        im2 = handles.images2(:,:,i);
        for j = 1:count
            sum1 = sum1 + im1(maskindex(j));
            sum2 = sum2 + im2(maskindex(j));
        end
        fprintf( fid, '%i,%.2f,%.2f\n',i,sum1,sum2 );
        fprintf('    %i   %.2f   %.2f\n',i,sum1,sum2 );
    end
    fclose(fid);
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dispimage1 = handles.images1(:,:,handles.nimages);
imagemax = max(dispimage1(:));
dispmask = dispimage1 * 0.0;
if ( ~isfield(handles,'disproi') )
    handles.disproi = createroimask2d(dispimage1);
else
    choice = questdlg('Redraw roi?','ROI choice','Yes','No','No');
    switch choice
        case 'Yes'
            handles.disproi = createroimask2d(dispimage1);
        case 'No'
    end
end

% Update handles structure
guidata(hObject, handles);
