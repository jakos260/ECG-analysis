function varargout = EGMplot(varargin)
% EGMPLOT M-file for EGMplot.fig
%      EGMPLOT, by itself, creates a new EGMPLOT or raises the existing
%      singleton*.
%
%      H = EGMPLOT returns the handle to a new EGMPLOT or the handle to
%      the existing singleton*.
%
%      EGMPLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EGMPLOT.M with the given input arguments.
%
%      EGMPLOT('Property','Value',...) creates a new EGMPLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EGMplot_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EGMplot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EGMplot

% Last Modified by GUIDE v2.5 24-Oct-2006 09:29:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EGMplot_OpeningFcn, ...
                   'gui_OutputFcn',  @EGMplot_OutputFcn, ...
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


% --- Executes just before EGMplot is made visible.
function EGMplot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EGMplot (see VARARGIN)

global SIG

pp=1;
sampF=1;
while pp<=length(varargin)
	if ischar(varargin{pp})
		key=lower(varargin{pp});
		switch key
			case 'sampt'
				sampT=varargin{pp+1};pp=pp+2;
			case 'sampf'
				sampT=1/varargin{pp+1};pp=pp+2;
			otherwise
				error('unknown parameter');
		end
	else
		if length(varargin)==1
			SIG.sig=varargin{1};
			SIG.nsig=1;
		elseif length(varargin)==2
			SIG.sig=varargin{2};
			SIG.sigt=varargin{1};
		else
			eval(['SIG.sig' num2str(pp) '=varargin{ ' num2str(pp) '};']);
			SIG.nsig=pp;
			SIG.sig=SIG.sig1;
		end
		pp=pp+1;
	end
end
if ~isfield(SIG,'t')
	SIG.sigt=(0:length(SIG.sig1)-1)*sampT;
end
SIG.Amin=floor(10*min(SIG.sig))/10;
SIG.Amax=ceil(10*max(SIG.sig))/10;
SIG.maxt=max(SIG.sigt);

SIG.zoom=0.1*SIG.maxt; if SIG.zoom > 1, SIG.zoom=round(SIG.zoom); end;
SIG.sigtmax=SIG.zoom;
SIG.sigtmin=min(SIG.sigt);
SIG.measure=0;
handles = guihandles(hObject);

set(hObject,'MenuBar','figure');
set(hObject,'ToolBar','figure');

% Choose default command line output for EGMplot
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EGMplot wait for user response (see UIRESUME)
% uiwait(handles.EGMplot);
set(handles.tmin_edit,'string',num2str(SIG.sigtmin));
set(handles.tzoom_edit,'string',num2str(SIG.zoom));
set(handles.Amin_edit,'string',num2str(SIG.Amin));
set(handles.Amax_edit,'string',num2str(SIG.Amax));
set(handles.checkbox_measure,'Value',SIG.measure);
set(handles.EGM_axes,'Ylim',[SIG.Amin SIG.Amax]);
zoomstep=SIG.zoom/SIG.maxt;
set(handles.xslider,'Sliderstep',[0.1*zoomstep,zoomstep],'Min',0,'Max',SIG.sigt(end)-zoomstep*SIG.sigt(end),'Value',0);
update_controls(handles);

% --- Outputs from this function are returned to the command line.
function varargout = EGMplot_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%--------------------------------------------------------------------------
function xslider_Callback(hObject, eventdata, handles)
% hObject    handle to xslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% global SIG

xmin=get(hObject,'Value');%*SIG.maxt 
if xmin > 1, xmin=round(xmin); end;
set(handles.tmin_edit,'String',num2str(xmin));
update_controls(handles);

%--------------------------------------------------------------------------

function analyse_axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to V_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global SIG


point1=get(hObject,'currentpoint');
delete(findobj('Tag','pline1')); delete(findobj('Tag','pline2')); delete(findobj('Tag','pline3'));
delete(findobj('Tag','textpl'));
line([point1(1,1) point1(1,1)],[SIG.Amin SIG.Amax],'Color','r','parent',hObject,'Tag','pline1','hittest','off');
finalRect = rbbox;                   % return figure units
point2 = get(hObject,'CurrentPoint');% button up detected
line([point2(1,1) point2(1,1)],[SIG.Amin SIG.Amax],'Color','r','parent',hObject,'Tag','pline2','hittest','off');
line([point1(1,1) point2(1,1)],[0 0],'Color','r','parent',hObject,'Tag','pline3','hittest','off');
if SIG.measure
	text(point1(1,1),SIG.Amax*0.9,[num2str(1000*(abs(point1(1,1)-point2(1,1))),4) 'ms'],'parent',hObject,'Tag','textpl','hittest','off')
elseif exist('AnalyseWindow')
	AnalyseWindow(SIG,point1(1,1),point2(1,1));
else
	xmax=max(point1(1,1),point2(1,1));
	xmin=min(point1(1,1),point2(1,1));
	ymax=max(point1(1,2),point2(1,2));
	ymin=min(point1(1,2),point2(1,2));	
	disp(num2str([xmin xmax ymin ymax],4));
	if xmax > xmin
		if xmax > SIG.sigt(end), xmax=SIG.sigt(end); xmin=xmax-SIG.zoom; 	set(handles.tmin_edit,'String',num2str(xmin));end
		SIG.zoom=max(xmax-xmin,1);
		set(handles.tzoom_edit,'String',num2str(SIG.zoom));
		set(handles.tmin_edit,'String',num2str(xmin));
		SIG.sigtmax=xmax;
		zoomstep=SIG.zoom/SIG.maxt;
		set(handles.xslider,'Sliderstep',[0.1*zoomstep,zoomstep]);
		set(handles.xslider,'Value',xmin);
		update_controls(handles);
	end
	if ymax > ymin
		SIG.Amin=ymin;
		SIG.Amax=ymax;
		set(handles.EGM_axes,'Ylim',[SIG.Amin SIG.Amax]);
		set(handles.Amin_edit,'string',num2str(SIG.Amin));
		set(handles.Amax_edit,'string',num2str(SIG.Amax));
	end
	if xmin==xmax && ymin==ymax
		zoom=SIG.maxt;

		SIG.zoom=zoom;
		SIG.sigtmin=0;
		set(handles.tmin_edit,'String','0');
		set(handles.tzoom_edit,'String',num2str(zoom,4));
		SIG.sigtmax=zoom;
		zoomstep=SIG.zoom/SIG.maxt;
		set(handles.xslider,'Sliderstep',[0.1*zoomstep,zoomstep]);
		set(handles.xslider,'Value',0);
		SIG.Amin=min(SIG.sig);
		SIG.Amax=max(SIG.sig);
		set(handles.EGM_axes,'Ylim',[SIG.Amin SIG.Amax]);
		set(handles.Amin_edit,'string',num2str(SIG.Amin));
		set(handles.Amax_edit,'string',num2str(SIG.Amax));

		update_controls(handles);
		
	end
end

%--------------------------------------------------------------------------
function tzoom_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tzoom_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global SIG

zoom=str2double(get(hObject,'String'));
if zoom > 0
	SIG.zoom=zoom;
	xmin=str2double(get(handles.tmin_edit,'String'));
	xmax=xmin+SIG.zoom;
	if xmax > SIG.sigt(end), xmax=SIG.sigt(end); xmin=xmax-SIG.zoom; 	set(handles.tmin_edit,'String',num2str(xmin));end
	SIG.sigtmax=xmax;
	zoomstep=SIG.zoom/SIG.maxt;
	set(handles.xslider,'Sliderstep',[0.1*zoomstep,zoomstep]);
	set(handles.xslider,'Value',xmin);
	update_controls(handles);
end

%--------------------------------------------------------------------------

function tmin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tmin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global SIG

xmin=str2double(get(hObject,'String'));
if xmin >=0
	if xmin < 0, set(hObject,'String','0');xmin=0; end;
	SIG.sigtmin=xmin;
	xmax=xmin+SIG.zoom;
	if xmax > SIG.sigt(end), xmax=SIG.sigt(end); xmin=xmax-SIG.zoom; 	set(handles.tmin_edit,'String',num2str(xmin));end
	SIG.sigtmax=xmax;
	zoomstep=SIG.zoom/SIG.maxt;
	set(handles.xslider,'Sliderstep',[0.1*zoomstep,zoomstep]);
	set(handles.xslider,'Value',xmin);
	update_controls(handles);
else
	set(handles.tmin_edit,'String',num2str(SIG.sigtmin));
end

%--------------------------------------------------------------------------
function Amax_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Amax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global SIG

SIG.Amax=str2double(get(hObject,'String'));
set(handles.EGM_axes,'Ylim',[SIG.Amin SIG.Amax])

%--------------------------------------------------------------------------
function Amin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Amin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global SIG

SIG.Amin=str2double(get(hObject,'String'));
set(handles.EGM_axes,'Ylim',[SIG.Amin SIG.Amax])

%==========================================================================

function update_controls(handles)

global SIG


posmin=round(str2double(get(handles.tmin_edit,'String')));
posmax=posmin+round(str2double(get(handles.tzoom_edit,'String')));
delete(get(handles.EGM_axes,'Children'));
select=find(SIG.sigt>=posmin & SIG.sigt <= posmax);
cols=['''b''';'''r''';'''k''';'''g''';'''c''';'''m''';'''y'''];
styles=['k- ';'r- ';'k--' ; 'r--' ;'k: ' ; 'k-.'];
Amin=floor(10*min(SIG.sig(select)))/10;
Amax=ceil(10*max(SIG.sig(select)))/10;
set(handles.EGM_axes,'Ylim',[Amin Amax])
if SIG.nsig>1
	for i=1:SIG.nsig
		eval(['line(SIG.sigt(select),SIG.sig' num2str(i) '(select),  ''Color'' , ' cols(i,:)  ', ''Parent'' ,handles.EGM_axes);']);
	end
else
	line(SIG.sigt(select),SIG.sig(select),'Color','b','Parent',handles.EGM_axes,'hittest','off');
end
set(handles.EGM_axes,'XLim',[posmin posmax]);
if SIG.zoom >= SIG.maxt
	set(handles.xslider,'Visible','off');
else
	set(handles.xslider,'Visible','on');
end
%==========================================================================

% --- Executes during object creation, after setting all properties.
function Amax_edit_CreateFcn(hObject, eventdata, handles)

if ispc,    
	set(hObject,'BackgroundColor','white');
else		   
	set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor')); 
end

% --- Executes during object creation, after setting all properties.
function Amin_edit_CreateFcn(hObject, eventdata, handles)

if ispc,    
	set(hObject,'BackgroundColor','white');
else		   
	set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor')); 
end

% --- Executes during object creation, after setting all properties.
function tmin_edit_CreateFcn(hObject, eventdata, handles)

if ispc,    
	set(hObject,'BackgroundColor','white');
else,		   
	set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor')); 
end

% --- Executes during object creation, after setting all properties.
function tzoom_edit_CreateFcn(hObject, eventdata, handles)

if ispc,    
	set(hObject,'BackgroundColor','white');
else		   
	set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor')); 
end

% --- Executes during object creation, after setting all properties.
function xslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





% --- Executes on button press in checkbox_measure.
function checkbox_measure_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_measure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_measure

global SIG

SIG.measure=get(hObject,'Value');


% --- Executes on button press in pushZoomReset.
function pushZoomReset_Callback(hObject, eventdata, handles)
% hObject    handle to pushZoomReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SIG

zoom=SIG.maxt;

SIG.zoom=zoom;
SIG.sigtmin=0;
set(handles.tmin_edit,'String','0');
set(handles.tzoom_edit,'String',num2str(zoom,4));
SIG.sigtmax=zoom;
zoomstep=SIG.zoom/SIG.maxt;
set(handles.xslider,'Sliderstep',[0.1*zoomstep,zoomstep]);
set(handles.xslider,'Value',0);
SIG.Amin=min(SIG.sig);
SIG.Amax=max(SIG.sig);
set(handles.EGM_axes,'Ylim',[SIG.Amin SIG.Amax]);
set(handles.Amin_edit,'string',num2str(SIG.Amin));
set(handles.Amax_edit,'string',num2str(SIG.Amax));
update_controls(handles);

