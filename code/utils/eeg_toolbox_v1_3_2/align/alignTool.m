function varargout = alignTool(varargin)
% ALIGNTOOL M-file for alignTool.fig
%      ALIGNTOOL, by itself, creates a new ALIGNTOOL or raises the existing
%      singleton*.
%
%      H = ALIGNTOOL returns the handle to a new ALIGNTOOL or the handle to
%      the existing singleton*.
%
%      ALIGNTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALIGNTOOL.M with the given input arguments.
%
%      ALIGNTOOL('Property','Value',...) creates a new ALIGNTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before alignTool_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to alignTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help alignTool

% Last Modified by GUIDE v2.5 26-Sep-2005 15:00:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @alignTool_OpeningFcn, ...
                   'gui_OutputFcn',  @alignTool_OutputFcn, ...
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


% --- Executes just before alignTool is made visible.
function alignTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to alignTool (see VARARGIN)

% Choose default command line output for alignTool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes alignTool wait for user response (see UIRESUME)
% uiwait(handles.figMain);


% --- Outputs from this function are returned to the command line.
function varargout = alignTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pbLoadFile.
function pbLoadFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoadFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

loadFile('',handles);


function loadFile(filename,handles)
%LOADFILE - Opens a dialog to open a file if a file is not provided.
%
%

% see if must open a dialog
if isempty(filename) 
  %[filename,pathname] = uigetfile({'*.*','All Files (*.*)'},'Pick a sync-pulse file','MultiSelect','on');
  [filename,pathname] = uigetfile('*','Pick a sync-pulse file','MultiSelect','on');
  %[filename,pathname] = uigetfile('*','Pick a sync-pulse file');
  if isequal(filename,0) | isequal(pathname,0)
    % canceled by user
    return
  end

  % temporary hack b/c of crappy matlab not having good standards
  if ischar(filename)
    filename = {filename};
  end
end

if length(filename) > 1
  % load two files
  filename1 = fullfile(pathname,filename{1});
  filename2 = fullfile(pathname,filename{2});
   
  % make sure the files exist
  if ~exist(filename1,'file')
    fprintf('ERROR: %s does not exist.\n',filename1);
  end
  if ~exist(filename2,'file')
    fprintf('ERROR: %s does not exist.\n',filename2);
  end
  
  % take diff of two files
  %dat = loadData(filename1) - loadData(filename2);
  [dat1,samplerate] = loadData(filename1);
  [dat2,samplerate] = loadData(filename2);
  dat = dat2-dat1;
  %dat = loadData(filename2) - loadData(filename1);
  
  % set the title
  set(handles.panelPulse,'Title',[filename1 ' - ' filename2]);
  
else
  % load one file
  filename = fullfile(pathname,filename{1});
  
  if ~exist(filename,'file')
    fprintf('ERROR: %s does not exist.\n',filename);
  end
 
  % load the file
  [dat,samplerate] = loadData(filename);
 
  % set the title
  set(handles.panelPulse,'Title',filename);
end


% get the axis data
pulseInfo = get(handles.axesMain,'UserData');

% plot it
cla
set(gca,'NextPlot','Add');
hDat = plot(handles.axesMain,dat);
hold on
hPulse = plot(handles.axesMain,[],'r*');
hold off

% create struct
pulseInfo.dat = dat;
pulseInfo.pulses = [];
pulseInfo.hDat = hDat;
pulseInfo.hPulse = hPulse;
pulseInfo.filename = filename;
pulseInfo.samplerate = samplerate;
  
% set the mode to Mark
pulseInfo = setMode('Mark',pulseInfo,handles);
set(handles.tbMark,'Value',1);
    
% save the struct
set(handles.axesMain,'UserData',pulseInfo);



function [dat,samplerate] = loadData(filename)
%LOADDATA - Load EEG data from a file
%

% get the data format
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(fileparts(filename));

% Open and load the file
fid = fopen(filename, 'r','l');
dat =  fread(fid, inf, dataformat);
fclose(fid);

% apply the gain
dat = dat.*gain;


%alignTool('axesButtonPress_Callback',gcbo,[],guidata(gcbo))


% --- Axes Button Press  ---
function axesButtonPress_Callback(hObject, eventdata, handles)
% Callback for button press in axes.
% Handles adding and removing of peaks.
%
%

% get the axis data
pulseInfo = get(handles.axesMain,'UserData');

minDistMS = 100;
%minDist = 10;
minDist = fix(minDistMS*pulseInfo.samplerate/1000);

% make sure not zooming or panning
if ~strcmp(pulseInfo.modeMarkPanZoom,'Mark') | isempty(pulseInfo)
  return
end

dat =  pulseInfo.dat;
axes(hObject);

% get the rectangle
point1 = get(gca,'CurrentPoint');    % button down detected
finalRect = rbbox;                   % return figure units
point2 = get(gca,'CurrentPoint');    % button up detected
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);

p = sortrows([point1;point2]);
p(p(:,1)<1,1) = 1;
p(p(:,1)>length(dat),1) = length(dat);


% set the xrange
xr = round(p(1,1)):round(p(2,1));

%get(handles.figMain,'SelectionType')

% see if adding or removing points
if strcmp(get(handles.figMain,'SelectionType'),'alt')
  % is shift, so delete
  % see if can delete anything
  if ~isempty(pulseInfo.hPulse)
    % get the current points
    xdat = get(pulseInfo.hPulse,'XData'); 
    ydat = get(pulseInfo.hPulse,'YData');
    
    % see if in box
    xind = find(ismember(xdat,xr));
    yind = find(ydat>=min(p(:,2)) & ydat<=max(p(:,2)));
    ind = intersect(xind,yind);
    
    % remove the points falling within the box
    xdat(ind) = [];
    ydat(ind) = [];
    
    % set the points back
    set(pulseInfo.hPulse,'XData',xdat);
    set(pulseInfo.hPulse,'YData',ydat);    
  end
else
  % adding points

  % get local max
  %[lmval,i] = lmax(dat(xr));
  i = findpeaks(dat(xr));
  
  if isempty(i)
    i = 1;
  end
  
  % set the global index points
  ind = min(xr) + i - 1;
  k = find(dat(ind)>=min(p(:,2)) & dat(ind)<=max(p(:,2)));
  ind = ind(k);
  
  if isempty(pulseInfo.hPulse)
    hold on
    pulseInfo.hPulse = plot(handles.axesMain,ind,dat(ind),'r*');
    hold off
  else
    % get the current points
    points = [get(pulseInfo.hPulse,'XData')' get(pulseInfo.hPulse,'YData')'];
    
    % append unique points
    ind = ind(find(~ismember(ind,points(:,1))));
    xdat = ind;
    ydat = dat(ind);
    points = [points ; xdat(:) ydat(:)];
    
    % sort by rows
    points = sortrows(points);
    
    % remove the second instance of any that are under min distance
    if minDist > 2
      mind = find(diff(points(:,1)) < minDist);
      if ~isempty(mind)
        points(mind+1,:) = [];
      end
    end

    
    % set the points back
    set(pulseInfo.hPulse,'XData',points(:,1),'YData',points(:,2));
  end
end

% save any modifications to the data
set(handles.axesMain,'UserData',pulseInfo);


% --- Executes on button press in tbZoom.
function tbZoom_Callback(hObject, eventdata, handles)
% hObject    handle to tbZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbZoom

% get the axis data
pulseInfo = get(handles.axesMain,'UserData');

% set the new mode
pulseInfo = setMode('Zoom',pulseInfo,handles);

% save any modifications to the data
set(handles.axesMain,'UserData',pulseInfo);


% --- Executes on button press in tbPan.
function tbPan_Callback(hObject, eventdata, handles)
% hObject    handle to tbPan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbPan

% get the axis data
pulseInfo = get(handles.axesMain,'UserData');

% set the new mode
pulseInfo = setMode('Pan',pulseInfo,handles);

% save any modifications to the data
set(handles.axesMain,'UserData',pulseInfo);


% --- Executes on button press in tbMark.
function tbMark_Callback(hObject, eventdata, handles)
% hObject    handle to tbMark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbMark

% get the axis data
pulseInfo = get(handles.axesMain,'UserData');

% set the new mode
pulseInfo = setMode('Mark',pulseInfo,handles);

% save any modifications to the data
set(handles.axesMain,'UserData',pulseInfo);


function pulseInfo = setMode(toggleMode, pulseInfo, handles)
% Handles the toggling of the three buttons.
% toggleMode can be one of {'Mark','Zoom','Pan'}
%
%

switch toggleMode
  case 'Mark'
    % set proper Zoom and Pan
    pan off;
    zoom off;
    
    % turn on the correct toggle
    set(handles.tbPan,'Value',0);
    set(handles.tbZoom,'Value',0);

  case 'Pan'
    % set proper Zoom and Pan
    pan on;
    zoom off;
    
    % turn on the correct toggle
    set(handles.tbMark,'Value',0);
    set(handles.tbZoom,'Value',0);
    
  case 'Zoom'
    % set proper Zoom and Pan
    pan off;
    zoom on;
    
    % turn on the correct toggle
    set(handles.tbMark,'Value',0);
    set(handles.tbPan,'Value',0);

end

pulseInfo.modeMarkPanZoom = toggleMode;



% --- Executes on button press in pbSavePulses.
function pbSavePulses_Callback(hObject, eventdata, handles)
% hObject    handle to pbSavePulses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% get the axis data
pulseInfo = get(handles.axesMain,'UserData');

if isempty(pulseInfo)
  return
end

% get the pulses
pulses = get(pulseInfo.hPulse,'XData');

% determine the suggested filename
if isstr(pulseInfo.filename)
  % is a single file, so just append sync.txt
  [fpath,fname,fext] = fileparts(pulseInfo.filename);
  outFileName = [fname fext '.sync.txt'];
else
  % is two files get base from both files to use
  [fpath1,fname1,fext1] = fileparts(pulseInfo.filename{1});
  [fpath2,fname2,fext2] = fileparts(pulseInfo.filename{2});
  outFileName = [fname1 fext1 fext2 '.sync.txt'];
end

if length(pulses) > 0
  % save to file
  [filename,pathname] = uiputfile(outFileName);
  if isequal(filename,0) | isequal(pathname,0)
    % canceled by user
    return
  end
  
  outfile = fullfile(pathname,filename);
  fid = fopen(outfile,'wt');
  fprintf(fid,'%d\n',pulses);
  fclose(fid);
  
end


% --- Executes on button press in pbFlipEEG.
function pbFlipEEG_Callback(hObject, eventdata, handles)
% hObject    handle to pbFlipEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% get the axis data
pulseInfo = get(handles.axesMain,'UserData');

% make sure there is data to flip
if isempty(pulseInfo)
  return
end

% flip the data
set(pulseInfo.hDat,'YData',-get(pulseInfo.hDat,'YData'));
set(pulseInfo.hPulse,'YData',-get(pulseInfo.hPulse,'YData'));
pulseInfo.dat = -pulseInfo.dat;
currYLim = get(handles.axesMain,'YLim');
set(handles.axesMain,'YLim',-[currYLim(2) currYLim(1)]);

% save any modifications to the data
set(handles.axesMain,'UserData',pulseInfo);


% --- Executes on selection change in lbEegFiles.
function lbEegFiles_Callback(hObject, eventdata, handles)
% hObject    handle to lbEegFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lbEegFiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lbEegFiles


% --- Executes during object creation, after setting all properties.
function lbEegFiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbEegFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lbBehFiles.
function lbBehFiles_Callback(hObject, eventdata, handles)
% hObject    handle to lbBehFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lbBehFiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lbBehFiles


% --- Executes during object creation, after setting all properties.
function lbBehFiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbBehFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tbPulse.
function tbPulse_Callback(hObject, eventdata, handles)
% hObject    handle to tbPulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbPulse

% make pulse visible and align not
set(handles.tbPulse,'Value',1);
set(handles.panelPulse,'Visible','On');

set(handles.tbAlign,'Value',0);
set(handles.panelAlign,'Visible','Off');


% --- Executes on button press in tbAlign.
function tbAlign_Callback(hObject, eventdata, handles)
% hObject    handle to tbAlign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbAlign

% make pulse not and align visible
set(handles.tbPulse,'Value',0);
set(handles.panelPulse,'Visible','Off');

set(handles.tbAlign,'Value',1);
set(handles.panelAlign,'Visible','On');


% --- Executes on button press in pbAddEEGFile.
function pbAddEEGFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbAddEEGFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

addFilesToListBox('Pick EEG pulse files',handles.lbEegFiles);


% --- Executes on button press in pbAddBehFile.
function pbAddBehFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbAddBehFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

addFilesToListBox('Pick behavioral pulse files',handles.lbBehFiles);



% --- Executes on button press in pbRemoveEegFile.
function pbRemoveEegFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbRemoveEegFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

removeFileFromListBox(handles.lbEegFiles);


% --- Executes on button press in pbRemoveBehFile.
function pbRemoveBehFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbRemoveBehFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

removeFileFromListBox(handles.lbBehFiles);



% --- Executes on button press in pbAddChanFile.
function pbAddChanFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbAddChanFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

addFilesToListBox('Pick sample channel files',handles.lbChanFiles);


% --- Executes on button press in pbRemoveChanFile.
function pbRemoveChanFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbRemoveChanFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

removeFileFromListBox(handles.lbChanFiles)




% --- Executes on button press in pbAddEventFile.
function pbAddEventFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbAddEventFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

addFilesToListBox('Pick Event files',handles.lbEventFiles);


% --- Executes on button press in pbRemoveEventFile.
function pbRemoveEventFile_Callback(hObject, eventdata, handles)
% hObject    handle to pbRemoveEventFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

removeFileFromListBox(handles.lbEventFiles);



%
function addFilesToListBox(dialogStr,lbHandle)
% Helper function 
%
%[filename,pathname] = uigetfile({'*.*','All Files (*.*)'},dialogStr,'MultiSelect','on');
[filename,pathname] = uigetfile('*',dialogStr,'MultiSelect','on');
%[filename,pathname] = uigetfile('*',dialogStr);
if isequal(filename,0) | isequal(pathname,0)
  % canceled by user
  return
end 

if ischar(filename)
  filename = {filename};
end

% add the new files
curFiles = get(lbHandle,'String');
if ~iscell(curFiles)
  curFiles = {};
end
ind = length(curFiles);
for f = 1:length(filename)
  % append the new file
  newfile = fullfile(pathname,filename{f});
  ind = ind + 1;
  curFiles{ind} = newfile;
end

% set it back with the new list
set(lbHandle,'String',curFiles);


%
function removeFileFromListBox(lbHandle)
% Helper function
%

% get the selected one
ind = get(lbHandle,'Value');
if ind > 0
  % remove that one
  curFiles = get(lbHandle,'String');
  if ~iscell(curFiles)
    curFiles = {};
    return
  end
  curFiles = curFiles([setdiff(1:length(curFiles),ind)]);
  
  % make sure the value will not be illegal
  if ind > length(curFiles)-1
    set(lbHandle,'Value',ind-1);
  end
  
  % set it back with the new list
  set(lbHandle,'String',curFiles);
end



% --- Executes on selection change in lbChanFiles.
function lbChanFiles_Callback(hObject, eventdata, handles)
% hObject    handle to lbChanFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lbChanFiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lbChanFiles


% --- Executes during object creation, after setting all properties.
function lbChanFiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbChanFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in cbIsFrei.
function cbIsFrei_Callback(hObject, eventdata, handles)
% hObject    handle to cbIsFrei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbIsFrei


% --- Executes on selection change in lbEventFiles.
function lbEventFiles_Callback(hObject, eventdata, handles)
% hObject    handle to lbEventFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lbEventFiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lbEventFiles


% --- Executes during object creation, after setting all properties.
function lbEventFiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbEventFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtTimeField_Callback(hObject, eventdata, handles)
% hObject    handle to txtTimeField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTimeField as text
%        str2double(get(hObject,'String')) returns contents of txtTimeField as a double


% --- Executes during object creation, after setting all properties.
function txtTimeField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTimeField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbRunAlign.
function pbRunAlign_Callback(hObject, eventdata, handles)
% hObject    handle to pbRunAlign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% gather all the input data
%srate_str = get(handles.popSamplerate,'String');
%ind = get(handles.popSamplerate,'Value');
%samplerate = str2num(srate_str{ind});

fprintf('\n\nRunning Alignment...\n\n');

beh_file = get(handles.lbBehFiles,'String');
eeg_file = get(handles.lbEegFiles,'String');
chan_files = get(handles.lbChanFiles,'String');
log_files = get(handles.lbEventFiles,'String');

% get the samplerate
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(fileparts(chan_files{1}));
fprintf('Samplerate: %g\n',samplerate);

ms_field = get(handles.txtTimeField,'String');

isfrei = get(handles.cbIsFrei,'Value');

% call the function
runAlign(samplerate,beh_file,eeg_file,chan_files,log_files,ms_field,isfrei)

fprintf('\nDone!!!\n\n');


% --- Executes on button press in pbEegFileMoveUp.
function pbEegFileMoveUp_Callback(hObject, eventdata, handles)
% hObject    handle to pbEegFileMoveUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

moveUpFileInListBox(handles.lbEegFiles);



% --- Executes on button press in pbEegFileMoveDown.
function pbEegFileMoveDown_Callback(hObject, eventdata, handles)
% hObject    handle to pbEegFileMoveDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


moveDownFileInListBox(handles.lbEegFiles);





% --- Executes on button press in pbChanFileMoveUp.
function pbChanFileMoveUp_Callback(hObject, eventdata, handles)
% hObject    handle to pbChanFileMoveUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

moveUpFileInListBox(handles.lbChanFiles);


% --- Executes on button press in pbChanFileMoveDown.
function pbChanFileMoveDown_Callback(hObject, eventdata, handles)
% hObject    handle to pbChanFileMoveDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

moveDownFileInListBox(handles.lbChanFiles);



%
function moveUpFileInListBox(lbHandle)
% Helper function
%

% get the selected one
ind = get(lbHandle,'Value');
if ind > 1
  % Flip with previous one
  curFiles = get(lbHandle,'String');
  if ~iscell(curFiles)
    curFiles = {};
    return
  end
  tempStr = curFiles{ind-1};
  curFiles{ind-1} = curFiles{ind};
  curFiles{ind} = tempStr;
 
  % set the new value
  set(lbHandle,'Value',ind-1);
  
  % set it back with the new list
  set(lbHandle,'String',curFiles);
end


%
function moveDownFileInListBox(lbHandle)
% Helper function
%

% get the selected one
ind = get(lbHandle,'Value');
curFiles = get(lbHandle,'String');
if ~iscell(curFiles)
  curFiles = {};
  return
end
if ind > 0 & length(curFiles) > ind
  % Flip with next one
  
  tempStr = curFiles{ind+1};
  curFiles{ind+1} = curFiles{ind};
  curFiles{ind} = tempStr;
 
  % set the new value
  set(lbHandle,'Value',ind+1);
  
  % set it back with the new list
  set(lbHandle,'String',curFiles);
end

