function varargout = selectBSMbeat(varargin)
% SELECTBSMBEAT M-file for selectBSMbeat.fig
%      SELECTBSMBEAT, by itself, creates a new SELECTBSMBEAT or raises the existing
%      singleton*.
%
%      H = SELECTBSMBEAT returns the handle to a new SELECTBSMBEAT or the handle to
%      the existing singleton*.
%
%      SELECTBSMBEAT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTBSMBEAT.M with the given input arguments.
%
%      SELECTBSMBEAT('Property','Value',...) creates a new SELECTBSMBEAT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before selectBSMbeat_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to selectBSMbeat_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help selectBSMbeat

% Last Modified by GUIDE v2.5 02-Jun-2011 09:55:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @selectBSMbeat_OpeningFcn, ...
                   'gui_OutputFcn',  @selectBSMbeat_OutputFcn, ...
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


%% --- Executes just before selectBSMbeat is made visible.
function selectBSMbeat_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to selectBSMbeat (see VARARGIN)

% global handles
set(handles.selectBSMFig,'Toolbar','figure')
set(handles.selectBSMFig,'menubar','figure')	

if length(varargin) < 1
	error('This routine needs at least two parameters');
else
	handles.DATA.ORG = varargin{1};
	handles.OPTIONS.funscal=2.5;
	handles.OPTIONS.zeromean=1;
	handles.OPTIONS.col='blue';
	handles.OPTIONS.linew=1;	
	handles.OPTIONS.label='';
	handles.DATA.sampT=1/1000;
	pp=2;
	while pp<=length(varargin)
		if ischar(varargin{pp})
			key=lower(varargin{pp});
			switch key
				case 'funscal'
					handles.OPTIONS.funscal=varargin{pp+1};pp=pp+2;
				case 'lay'
					handles.DATA.LAY=varargin{pp+1};pp=pp+2;
				case 'remove'
					handles.DATA.remove=varargin{pp+1};pp=pp+2;
				case 'sampt'
					handles.DATA.sampT=varargin{pp+1};pp=pp+2;
				case 'color'
					handles.OPTIONS.col=varargin{pp+1};pp=pp+2;
				case 'linewidth'
					handles.OPTIONS.linew=varargin{pp+1};pp=pp+2;
				case 'label'
					handles.OPTIONS.label=varargin{pp+1};pp=pp+2;
					
				otherwise
					error('unknown parameter');
			end
		end
	end
end
sumData = sum(handles.DATA.ORG,2);
handles.DATA.remove = sumData > mean(sumData) + std(sumData);
handles.DATA.PQint = 180; % default PQinterval + a part of the QRS interval
handles.DATA.len20ms=ceil(0.020./handles.DATA.sampT);

handles.OPTIONS.viewt0=1;
handles.OPTIONS.viewt1=min(handles.DATA.len20ms + 4/(handles.DATA.sampT),size(handles.DATA.ORG,2));

% Choose default command line output for selectBSMbeat
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes selectBSMbeat wait for user response (see UIRESUME)
% uiwait(handles.selectBSMFig);

set(handles.sliderSignal,'sliderStep',[0.01,0.1])
set(handles.sliderSignal,'min',1);
set(handles.sliderSignal,'max',size(handles.DATA.ORG,2))
set(handles.sliderSignal,'value',1);
set(handles.checkboxZeromean,'value',handles.OPTIONS.zeromean);

% set(handles.RMSaxes,'ButtonDownFcn',['selectBSMbeat(''RMSaxes_ButtonDownFcn'' ,gcbo,[],guidata(gcbo))']);
% set(handles.plotAxes,'ButtonDownFcn','selectBSMbeat(''plotAxes_ButtonDownFcn'' ,gcbo,[],guidata(gcbo))');

doBaseline(handles);

%% --- Outputs from this function are returned to the command line.
function varargout = selectBSMbeat_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
if nargout > 1
    varargout{2} = handles.DATA;
end

%%========================================================================
function doBaseline(handles)

handles.DATA.BSM = handles.DATA.ORG;

handles.DATA.BSM(handles.DATA.remove==1,:) = 0;
% lowpass at 200 Hz moving average
handles.DATA.BSM_low200 = lowpassma(handles.DATA.BSM,1/(200*handles.DATA.sampT));

% bandpass at 0.1 to 40 Hz moving average
BSM = lowpassma(handles.DATA.BSM,1/(40*handles.DATA.sampT)); 
handles.DATA.BSMfilt = BSM - ( lowpassma(BSM, 1 / ( 0.1 * handles.DATA.sampT ) ) );

% do zeromean
if handles.OPTIONS.zeromean==1
    handles.DATA.BSMfilt = zeromean(handles.DATA.BSMfilt);    
    handles.DATA.BSM_low200 = zeromean(handles.DATA.BSM_low200);    
end
handles.DATA.BSMOUT = zeromean(handles.DATA.BSMfilt); % OUTPUT BY DEFINITION ZEROMEAN




%% find baseline points --------------------------------------------------
rrms = rms(handles.DATA.BSMfilt);
rrms = diffrows(rrms);

%     plot(rrms)
qrss=find(rrms>0.5*max(rrms));
qrss(qrss < 300 )=[];

handles.DATA.BEATS = qrss(1) - handles.DATA.PQint;
for i=2:length(qrss)
    if qrss(i) - qrss(i-1) > 400 % minimal interval between qrss
        handles.DATA.BEATS = [handles.DATA.BEATS; qrss(i) - handles.DATA.PQint];
    end
end

if ~isempty(handles.DATA.BEATS)
    zerotims=zeros(length(handles.DATA.BEATS)+2,1);
    zerotims(1) = 1;
    zerotims(2:end-1) = handles.DATA.BEATS;
    zerotims(end) = length(rrms);
    BASEL=spline(zerotims,handles.DATA.BSM_low200(:,zerotims),1:size(handles.DATA.BSM_low200,2));
    handles.DATA.BSM_low200 = handles.DATA.BSM_low200 - BASEL;   
end

handles.DATA.RMS     = rms(handles.DATA.BSM_low200);
handles.DATA.RMSfilt = rms(handles.DATA.BSMfilt);

if ~isfield(handles.DATA,'SELBEATS') || ...
    length(handles.DATA.SELBEATS) ~= length(handles.DATA.BEATS)
    % resize SELBEATS, i.e. structure is cleared
    handles.DATA.SELBEATS = zeros(size(handles.DATA.BEATS));
end
% Update handles structure
guidata(handles.selectBSMFig, handles);


plotBSM(handles);

%%========================================================================
function plotBSM(handles)

% do highpass
if get(handles.checkboxHighPass,'Value') == 1
    bsm = handles.DATA.BSMfilt;
    rrms = handles.DATA.RMSfilt;
else
    bsm = handles.DATA.BSM_low200;
    rrms = handles.DATA.RMS;
end

set(handles.axesBeats,'Visible','on');
set(handles.selectBSMFig,'CurrentAxes',handles.axesBeats);
delete(get(handles.axesBeats,'Children'));	   

rtop0=round(300/(1000*handles.DATA.sampT));%find(rms(beat0)==max(rms(beat0)));
nr=0;
for k=1:length(handles.DATA.BEATS)-1
    if handles.DATA.BEATS(k) >= handles.OPTIONS.viewt0 && handles.DATA.BEATS(k+1) <= handles.OPTIONS.viewt1
        nr=nr+1;
        beat = bsm(:,handles.DATA.BEATS(k):handles.DATA.BEATS(k+1));
        rtop = find(rms(beat)==max(rms(beat)));
        if rtop<rtop0
            beat=[zeros(size(beat,1),abs(rtop-rtop0)+1) beat];
        elseif rtop>rtop0
            beat=beat(:,abs(rtop-rtop0):end);
        end
        plot(rms(beat),'r','parent',handles.axesBeats);
        hold on
    end
    axis tight
end
axis off

%% bsm signals
set(handles.selectBSMFig,'CurrentAxes',handles.plotAxes);
delete(get(handles.plotAxes,'Children'));

nplts=size(handles.DATA.LAY,1) - 1; % note: nplts may be smaller than nlds
raster=handles.DATA.LAY(1,:);
axis([0 raster(1) 0 2*raster(2)] )
axis manual
% axis('off')
box off
funscal=handles.OPTIONS.funscal;
hold on

nrastercols=handles.DATA.LAY(1,1);

bsm  = bsm(:,max(1,handles.OPTIONS.viewt0):min(size(bsm,2),handles.OPTIONS.viewt1));
t = 0 : size(bsm,2) - 1;
t = t*.9/max(t) * nrastercols/(nrastercols);

for i=2:nplts+1
	j=handles.DATA.LAY(i,3);
	yshift(j) = 2*raster(2)+1 - 2*handles.DATA.LAY(i,2);
	xshift(j) = handles.DATA.LAY(i,1)-1;
	if handles.DATA.remove(j)
		plot(t(ceil(length(t)/2))+xshift(j),yshift(j),'xr','Markersize',40,'hittest','off');
	else
		plot(t+xshift(j),funscal*bsm(j,:)+yshift(j),handles.OPTIONS.col,'linewidth',handles.OPTIONS.linew,'hittest','off')
		plot([t(1)+xshift(j) t(end)+xshift(j)],[yshift(j) yshift(j)],':k','hittest','off')
	end
	text(xshift(j)+.4,yshift(j)+.25,num2str(j));
end

extr = extremes(bsm);
plot(t(extr(3))+xshift(extr(2)),handles.OPTIONS.funscal*extr(1)+yshift(extr(2)),'*b');
plot(t(extr(6))+xshift(extr(5)),handles.OPTIONS.funscal*extr(4)+yshift(extr(5)),'*r');

xl=[0  1];
yl=[0  0];
set(line(xl,yl),'color','k')
text(0.05, 0.6,sprintf('%0.3f mV',1/handles.OPTIONS.funscal),'Tag','mVscale');
text(-handles.DATA.LAY(1,1)/10, -.25,sprintf('%d %s %0.3f %s',size(bsm,2),' samples; ',(size(bsm,2) - 1)*handles.DATA.sampT,'s'),'Tag','timetext');

%% RMS 
set(handles.selectBSMFig,'CurrentAxes',handles.RMSaxes);
delete(get(handles.RMSaxes,'Children'));

maxA=max(rrms);
minA=min(rrms);
rrms = rrms(:,max(1,handles.OPTIONS.viewt0):min(size(bsm,2),handles.OPTIONS.viewt1));
rmst = 0 : length(rrms)-1;

plot(rmst,rrms,'k','hittest','off','parent',handles.RMSaxes);

hold on
beats = handles.DATA.BEATS;
beats = beats - handles.OPTIONS.viewt0;
beats(beats<0)=[];
beats(beats > length(rrms))=[];
for i=1:length(handles.DATA.SELBEATS)
    if handles.DATA.SELBEATS(i)==1 && ...
       handles.DATA.BEATS(i) > handles.OPTIONS.viewt0 && handles.DATA.BEATS(i+1) < handles.OPTIONS.viewt1
        plot(handles.DATA.t(handles.DATA.BEATS(i):handles.DATA.BEATS(i+1)),zeros(size(handles.DATA.t(handles.DATA.BEATS(i):handles.DATA.BEATS(i+1)))),'r','linewidth',10,'parent',handles.RMSaxes);
        hold on;
    end
end

plot(rmst(beats),rrms(beats),'or','Markersize',10,'parent',handles.RMSaxes);
axis([rmst(1) rmst(end) minA maxA]);

% set(handles.RMSaxes ,'ButtonDownFcn','selectBSMbeat(''RMSaxes_ButtonDownFcn'' ,gcbo,[],guidata(gcbo))');
% set(handles.plotAxes,'ButtonDownFcn','selectBSMbeat(''plotAxes_ButtonDownFcn'' ,gcbo,[],guidata(gcbo))');


%% --- Executes on slider movement.
function sliderSignal_Callback(hObject, eventdata, handles)
% hObject    handle to sliderSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


val=round(get(handles.sliderSignal,'Value'));
viewwin=round(handles.OPTIONS.viewt1-handles.OPTIONS.viewt0);
if viewwin+val>size(handles.DATA.BSM,2)-handles.DATA.len20ms
	handles.OPTIONS.viewt1=size(handles.DATA.BSM,2)-handles.DATA.len20ms;
	handles.OPTIONS.viewt0=handles.OPTIONS.viewt1-viewwin;
else
	handles.OPTIONS.viewt0=round(val);
	handles.OPTIONS.viewt1=round(val+viewwin);
end
% Update handles structure
guidata(hObject, handles);

plotBSM(handles)



%% --- Executes when selectBSMFig is resized.
function selectBSMFig_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to selectBSMFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles)
	return
end
figp=get(handles.selectBSMFig,'Position');
plotp=get(handles.plotAxes,'Position');
rmsp=get(handles.RMSaxes,'Position');
contp=get(handles.uipanelControls,'Position');
contp(1)=.5;contp(2)=.5;
y=contp(4)+5;
dy=figp(4)-(contp(2)+y);

rmsp(1)=0.5;rmsp(2)=y;rmsp(3)=figp(3)-1;rmsp(4)=0.349*dy;
plotp(1)=0.5;plotp(2)=0.35*dy+y;plotp(3)=figp(3)-1;plotp(4)=0.65*dy;

set(handles.plotAxes,'Position',plotp);
set(handles.RMSaxes,'Position',rmsp);
set(handles.uipanelControls,'Position',contp);


%% --- Executes on button press in pushZoomIn.
function pushZoomIn_Callback(hObject, eventdata, handles)
% hObject    handle to pushZoomIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val=round((handles.OPTIONS.viewt1-handles.OPTIONS.viewt0)/2);
viewwin=round((handles.OPTIONS.viewt1-handles.OPTIONS.viewt0)/2);
if viewwin/2+val>size(handles.DATA.BSM,2)-handles.DATA.len20ms
	handles.OPTIONS.viewt1=size(handles.DATA.BSM,2)-handles.DATA.len20ms;
	handles.OPTIONS.viewt0=handles.OPTIONS.viewt1-viewwin;
else
	handles.OPTIONS.viewt0=val-round(viewwin/2);
	handles.OPTIONS.viewt1=val+round(viewwin/2);
end
updateSignalSlider(handles)
% Update handles structure
guidata(hObject, handles);

plotBSM(handles)


%% --- Executes on button press in pushZoomOut.
function pushZoomOut_Callback(hObject, eventdata, handles)
% hObject    handle to pushZoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




val=round((handles.OPTIONS.viewt1-handles.OPTIONS.viewt0)/2);
viewwin=(handles.OPTIONS.viewt1-handles.OPTIONS.viewt0)*2;
if viewwin/2+val>size(handles.DATA.BSM,2)-handles.DATA.len20ms
	handles.OPTIONS.viewt1=size(handles.DATA.BSM,2)-handles.DATA.len20ms;
	handles.OPTIONS.viewt0=max(1,handles.OPTIONS.viewt1-viewwin);
elseif val-viewwin/2<handles.DATA.len20ms
	handles.OPTIONS.viewt0=handles.DATA.len20ms;
	handles.OPTIONS.viewt1=val+viewwin;
else
	handles.OPTIONS.viewt0=val-viewwin/2;
	handles.OPTIONS.viewt1=val+viewwin/2;
end
updateSignalSlider(handles)
% Update handles structure
guidata(hObject, handles);

plotBSM(handles)


function updateSignalSlider(handles)

slidstep=(handles.OPTIONS.viewt1-handles.OPTIONS.viewt0)/size(handles.DATA.BSM,2);
val= get(handles.sliderSignal,'value')
set(handles.sliderSignal,'Sliderstep',[0.1*slidstep slidstep],'Value',handles.OPTIONS.viewt0);
set(handles.sliderSignal,'value',val)
if handles.OPTIONS.viewt0==1 && handles.OPTIONS.viewt1 == size(handles.DATA.BSM,2)
    set(handles.sliderSignal,'Enable','off' );
else
    set(handles.sliderSignal,'Enable','on' );
end


%% --- Executes on mouse press over axes background.
function RMSaxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to RMSaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

point1=get(hObject,'currentpoint');
selTime=round(point1(1)/handles.DATA.sampT);% + handles.OPTIONS.viewt0);
selMarker = find(abs(handles.DATA.BEATS-selTime) < 20);
maxA=max(handles.DATA.RMS);

if strcmp(get(handles.selectBSMFig,'SelectionType'),'alt')% add zeropoint
    
    markers = handles.DATA.BEATS-selTime;
    markers(markers<0) =10000;
    selBeat = find(markers == min(markers))-1;
    if selBeat > 1
        handles.DATA.SELBEATS(selBeat) = 1;
    end
    
	
elseif ~isempty(selMarker)  
    figpoint1=get(handles.selectBSMFig,'currentpoint');
    rect=[figpoint1(1,1), figpoint1(1,2), 2, maxA];   
    rbbox;
    point1=get(hObject,'currentpoint');
    handles.DATA.BEATS(selMarker) = round(point1(1)/handles.DATA.sampT)';% + handles.OPTIONS.viewt0;       
elseif 0

    delete(findobj('Tag','pline1')); delete(findobj('Tag','pline2')); 
    line([point1(1,1) point1(1,1)],[0 maxA],'Color','r','parent',hObject,'Tag','pline1','hittest','off');
    rbbox; 

    % return figure units
    point2 = get(hObject,'CurrentPoint');% button up detected

    line([point2(1,1) point2(1,1)],[0 maxA],'Color','r','parent',hObject,'Tag','pline2','hittest','off');
    xmax=max(point1(1,1),point2(1,1));
    xmin=min(point1(1,1),point2(1,1));
    if xmax < size(handles.DATA.BSM,2) - 1 %add one second min width
        xmax = max(xmin+1,xmax);
    else
        xmin = min(xmin,xmax-1);
    end
%     disp(num2str([xmin xmax],4));
    handles.OPTIONS.viewt0 = max(1,round(xmin/handles.DATA.sampT));
    handles.OPTIONS.viewt1 = min(size(handles.DATA.BSM,2),round(xmax/handles.DATA.sampT));
    delete(findobj('Tag','pline1')); delete(findobj('Tag','pline2'));     
    updateSignalSlider(handles);
%     set(handles.sliderSignal,'Value',(handles.OPTIONS.viewt0+handles.OPTIONS.viewt1)/2);
end
% Update handles structure
guidata(hObject, handles);

plotBSM(handles)



% --- Executes on mouse press over axes background.
function plotAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to plotAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% 	yshift(j)=2*raster(2)+1-2*handles.DATA.LAY(i,2);
% 	xshift(j)=handles.DATA.LAY(i,1)-1;

point1=get(hObject,'currentpoint');
x=floor(1+point1(1,1));
y=-round((point1(1,2)-handles.DATA.LAY(1,2)*2-1)/2);
% disp(num2str([point1(1,1) point1(1,2) x y]))
a=find(handles.DATA.LAY(:,1)==x & handles.DATA.LAY(:,2)==y);
if ~isempty(a)
	handles.DATA.remove(handles.DATA.LAY(a,3))=xor(handles.DATA.remove(handles.DATA.LAY(a,3)),1);
end
% Update handles structure
guidata(hObject, handles);

doBaseline(handles);

% --- Executes on button press in checkboxHighPass.
function checkboxHighPass_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxHighPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxHighPass


plotBSM(handles);



% --- Executes on button press in checkboxZeromean.
function checkboxZeromean_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxZeromean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



handles.OPTIONS.zeromean=get(hObject,'value');
% Update handles structure
guidata(hObject, handles);

doBaseline(handles);


% --- Executes on scroll wheel click while the figure is in focus.
function selectBSMFig_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to selectBSMFig (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)

handles.OPTIONS.funscal = max(0.5,handles.OPTIONS.funscal - eventdata.VerticalScrollCount/10.0);
% Update handles structure
guidata(hObject, handles);

plotBSM(handles);

% --- Executes on button press in recompBaseline_Button.
function recompBaseline_Button_Callback(hObject, eventdata, handles)
% hObject    handle to recompBaseline_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

doBaseline(handles);

% --- Executes when user attempts to close selectBSMFig.
function selectBSMFig_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to selectBSMFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
varargout{2} = handles.DATA;

delete(hObject);


% --- Executes during object creation, after setting all properties.
function selectBSMFig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectBSMFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called





% --- Executes during object creation, after setting all properties.
function sliderSignal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
