function varargout = torsopot(varargin)
% TORSOPOT M-file for torsopot.fig
%      TORSOPOT, by itself, creates a new TORSOPOT or raises the existing
%      singleton*.
%
%      H = TORSOPOT returns the handle to a new TORSOPOT or the handle to
%      the existing singleton*.
%
%      TORSOPOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TORSOPOT.M with the given input arguments.
%
%      TORSOPOT('Property','Value',...) creates a new TORSOPOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before torsopot_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to torsopot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help torsopot

% Last Modified by GUIDE v2.5 17-Jan-2008 10:23:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @torsopot_OpeningFcn, ...
                   'gui_OutputFcn',  @torsopot_OutputFcn, ...
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


% --- Executes just before torsopot is made visible.
function torsopot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to torsopot (see VARARGIN)

% Choose default command line output for torsopot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes torsopot wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global MYDAT

if length(varargin) < 3
	error('This routine needs at least 3 parameters');
else
	set(handles.figure1,'Toolbar','figure')
	set(handles.figure1,'menubar','figure')	
	
	MYDAT.minmax=[];
	MYDAT.ITRI=[];
	MYDAT.VER=[];
	MYDAT.DATA=[];
	MYDAT.ref=0;
	MYDAT.t0=1;

	pp=1;
	
	while pp<=length(varargin)
		if ischar(varargin{pp})
			key=lower(varargin{pp});
			switch key
				case 'faces'
					MYDAT.ITRI=varargin{pp+1};pp=pp+2;
				case 'vertices'
					MYDAT.VER=varargin{pp+1};pp=pp+2;
				case 'data'
					MYDAT.DATA=varargin{pp+1};pp=pp+2;
				case 'range'
					MYDAT.minmax=varargin{pp+1};pp=pp+2;
				case 'ref'
					MYDAT.ref=varargin{pp+1};pp=pp+2;
				case 'names'
					MYDAT.names=varargin{pp+1};pp=pp+2;
				otherwise
					error('unknown parameter');
			end
		end
	end
end
if isempty(MYDAT.VER) || isempty(MYDAT.ITRI) || isempty(MYDAT.DATA)
	error('missing necessary information to plot');
end

if length(size(MYDAT.DATA))
	MYDAT.DATAALL=MYDAT.DATA;
	MYDAT.DATA=squeeze(MYDAT.DATA(1,:,:));
	MYDAT.DATAALLindex=1;
end

[MYDAT.ADJ,MYDAT.DIST]=graphdist(MYDAT.ITRI,MYDAT.VER,4);
if isempty(MYDAT.minmax)
	MYDAT.minmax(1)=-max(max(max(abs(MYDAT.DATA))));
	MYDAT.minmax(2)= max(max(max(abs(MYDAT.DATA))));
end
if isfield(MYDAT,'names')
	set(handles.ListMaps,'String',MYDAT.names)
end
set(handles.datslider,'Sliderstep',[0.01 0.1],'Min',1,'Max',size(MYDAT.DATA,2),'Value',1);
colormap(loadmat('mymap_basic.mcm'))
update(handles)




%%
% --- Outputs from this function are returned to the command line.
function varargout = torsopot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function datslider_Callback(hObject, eventdata, handles)
% hObject    handle to datslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global MYDAT
MYDAT.t0=min(max(round(get(hObject,'Value')),1),size(MYDAT.DATA,2));
update(handles)

% --- Executes during object creation, after setting all properties.
function datslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%%
function selectnode(hObject, eventdata, handles)

global MYDAT

L   =   get(gca,'currentpoint');
ver=get(findobj(gca,'Type','patch'),'Vertices');
itri=get(findobj(gca,'Type','patch'),'faces');
if iscell(ver)
	ver=ver{end};
	itri=itri{end};
end
% intersections:
SELECT=linetris(ver,itri,L(1,:),L(2,:));
% only the facing ones are relavant
% SELECT=SELECT(SELECT(:,2)>0,:);

[small,is]=min(SELECT(:,5));
veris=itri(SELECT(is,1),:);
[small,ii]=min(sum((ver(veris,:)-(ones(1,3)'*(L(1,:)+SELECT(is,5)*(L(2,:)-L(1,:))))).^2'));
node=veris(ii);
delete(findobj('Tag','dot'));
line(ver(node,1),ver(node,2),ver(node,3),'Markersize',8,'Marker','o','Color','w','Linewidth',3,'Tag','dot');
% triag = SELECT(is,1);

MYDAT.ref=node;
update(handles)

%%
function update(handles)

global MYDAT

cla(handles.torsoaxis)
if MYDAT.ref~=0
	A=MYDAT.DATA(:,MYDAT.t0)-MYDAT.DATA(MYDAT.ref,MYDAT.t0);
	line(MYDAT.VER(MYDAT.ref,1),MYDAT.VER(MYDAT.ref,2),MYDAT.VER(MYDAT.ref,3),'Markersize',8,'Marker','o','Color','w','Linewidth',3,'Tag','dot');
	
else
	A=MYDAT.DATA(:,MYDAT.t0);
end
set(handles.figure1,'CurrentAxes',handles.torsoaxis);
MYDAT.hp=patch('Faces',MYDAT.ITRI,'Vertices',MYDAT.VER,...
	  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
	  'FaceVertexCData',A,'FaceColor','interp','Parent',handles.torsoaxis,...
	  'edgecolor','none','FaceAlpha',1,'buttondownFcn','torsopot(''selectnode'',gcbo,[],guidata(gcbo))');
view(90,0);
axis off equal tight
caxis([MYDAT.minmax(1) MYDAT.minmax(2)]);
contourlines(MYDAT.VER,MYDAT.ITRI,A,'sym',1,'extremes',1,'delta',0.2); 
colorbar
if MYDAT.ref
	distline(MYDAT.VER,MYDAT.ITRI,MYDAT.ref,45,MYDAT.ADJ)
end


cla(handles.ECGaxis)
plot(MYDAT.DATA(18,:),'parent',handles.ECGaxis,'Tag','myECG');
line([MYDAT.t0 MYDAT.t0],[-1 1],'color','r','parent',handles.ECGaxis,'Tag','hairline');

%% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global MYDAT
MYDAT.ref=0;
update(handles)


% --- Executes on selection change in ListMaps.
function ListMaps_Callback(hObject, eventdata, handles)
% hObject    handle to ListMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ListMaps contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListMaps

global MYDAT
a=get(hObject,'Value');
MYDAT.DATAALLindex=get(hObject,'Value');
MYDAT.DATA=squeeze(MYDAT.DATAALL(MYDAT.DATAALLindex,:,:));
rrms=rms(MYDAT.DATA);
MYDAT.t0=find(max(rrms)==rrms);
update(handles);






% --- Executes during object creation, after setting all properties.
function ListMaps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ListMaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


