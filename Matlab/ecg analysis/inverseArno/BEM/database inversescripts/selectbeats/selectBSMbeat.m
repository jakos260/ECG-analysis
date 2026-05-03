function varargout = selectBSMbeat(varargin)
% SELECTBSMBEAT M-file for selectBSMbeat.fig
%      SELECTBSMBEAT, by itself, creates a new SELECTBSMBEAT or raises the existing
%      singleton*.
%
%      H = SELECTBSMBEAT returns the handle to a new SELECTBSMBEAT or the handle to
%      the existing singleton*.DATA.BSM_low200
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

% oostep1:  DATA.PQint set to 120 ms instead of 180 samples.
%           round lowpassma window
%           Add ui control for DATA.PQint
%           PQint init at 220 ms, at 280 ms
% 20120622 oostep1: select channels for zeromean, ShowBSMonThorax, doBaseline when changing PQInt
% 20120813 oostep1: zoom to center.



% Last Modified by GUIDE v2.5 26-Jun-2013 15:49:25

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

global OPTIONS
global DATA

drawnow

% global handles
set(handles.selectBSMFig,'Toolbar','figure')
set(handles.selectBSMFig,'menubar','figure')

OPTIONS.funscal=0.5;
OPTIONS.zeromean=0;
OPTIONS.col='blue';
OPTIONS.linew=1;
OPTIONS.label='';
OPTIONS.leadsys='pigbsm';%'pigbsm'/'ams'; % default leadsystem for ShowBSMonThorax
OPTIONS.zmchannels=[];

if length(varargin) < 1
    % 	error('This routine needs at least two parameters');
    readmatbutton_Callback(handles.readmatbutton,eventdata,handles); % read DATA from file
else
    set(handles.Exportbutton,'Visible','off'); % export handled by caller
    DATA.ORG = varargin{1};
    DATA.sampT=1/1000;
    
%     % Add 3s of 0 if you need to select complexes at he beginning.
%     warning(' Added 3s of 0 to the beginning');
%     DATA.ORG=[zeros(size(DATA.ORG,1),3000) DATA.ORG];
%     DATA.RMS=[zeros(size(DATA.RMS,1),3000) DATA.RMS];
%     DATA.RMSfilt=[zeros(size(DATA.RMSfilt,1),3000) DATA.RMSfilt];
    
    pp=2;
    while pp<=length(varargin)
        if ischar(varargin{pp})
            key=lower(varargin{pp});
            switch key
                case 'funscal'
                    OPTIONS.funscal=varargin{pp+1};pp=pp+2;
                case 'lay'
                    DATA.LAY=varargin{pp+1};pp=pp+2;
                case 'remove'
                    DATA.remove=varargin{pp+1};pp=pp+2;
                case 'sampt'
                    DATA.sampT=varargin{pp+1};pp=pp+2;
                case 'color'
                    OPTIONS.col=varargin{pp+1};pp=pp+2;
                case 'linewidth'
                    OPTIONS.linew=varargin{pp+1};pp=pp+2;
                case 'label'
                    OPTIONS.label=varargin{pp+1};pp=pp+2;
                case 'leadsys'
                    if isempty(OPTIONS.leadsys)
                        pp=pp+2; % keep at default
                    else
                        OPTIONS.leadsys=varargin{pp+1};pp=pp+2;
                    end
                    
                otherwise
                    error('unknown parameter');
            end
        end
    end
    if isfield(DATA,'SELBEATS')
        DATA.SELBEATS =[];
    end   
    
    % Split off EGM channels, not in lay
    if size(DATA.ORG,1)>254 % Bucket hack
        DATA.EGM=DATA.ORG([1:31 96:size(DATA.ORG,1)],:);
        DATA.ORG=DATA.ORG(33:96,:);
    else
        maxlay=max(DATA.LAY(:,3));
        if maxlay<size(DATA.ORG,1)
            DATA.EGM=DATA.ORG(maxlay+1:end,:);
            DATA.ORG=DATA.ORG(1:maxlay,:);
        end
    end
    
    %Init OverlaySignal list. Do this before stripping EGM channels
    siglist=[{'(none)'};cellstr(num2str((1:size(DATA.ORG,1))'))];
    if isfield(DATA,'EGM');
        for i=1:size(DATA.EGM,1)
            siglist{end+1}=sprintf('egm%03d',i);
        end
    end
    set(handles.OverlaySignal,'String',siglist)
    set(handles.OverlaySignal,'Value',1);

    % BSM = lowpassma(DATA.ORG,round(1/(40*DATA.sampT)));
    sumData = mean(abs(DATA.ORG -  lowpassma(DATA.ORG, round( 1/(0.1*DATA.sampT)) ) ),2);
    % sumData = mean((DATA.ORG),2);
    DATA.remove = sumData > 7 * std(sumData);
    %     DATA.remove(33)=1 % temp for 059
%     DATA.remove([1 5 63 64])=1; % Pig10

    pqint=260; % PQ interval in ms
    
%     pqint=30 % PQ voor Sync Pig10
%     DATA.remove([1,5,25,33,39,45,57,63:68])=1 % Voor Pig09
    
    DATA.PQint = round(pqint/1000/DATA.sampT); % default PQinterval + a part of the QRS interval, in samples
    set(handles.PQint,'String',num2str(pqint));
    
    OPTIONS.viewt0=1;
    %     OPTIONS.viewt1=ceil(max(100,size(DATA.ORG,2)/4)); % PO
    OPTIONS.viewt1=ceil(min(5000,size(DATA.ORG,2)/2));
    
    %     OPTIONS.viewt0=size(DATA.ORG,2)/4 % hack for ajm059
    %     OPTIONS.viewt1=size(DATA.ORG,2)*3/4% hack for ajm059
    
    % Choose default command line output for selectBSMbeat
    handles.output = hObject;
    
    % Update handles structure
    guidata(hObject, handles);
    
    % UIWAIT makes selectBSMbeat wait for user response (see UIRESUME)
    % uiwait(handles.selectBSMFig);
    
    set(handles.sliderSignal,'sliderStep',[0.01,0.1])
    set(handles.sliderSignal,'min',1);
    
    %     set(handles.sliderSignal,'value',1);
    set(handles.sliderSignal,'value',size(DATA.ORG,2)/2); % hack for ajm059
    
    slidstep = ((OPTIONS.viewt1 - OPTIONS.viewt0) / size(DATA.ORG,2));
    set(handles.sliderSignal,'Sliderstep',[0.1*slidstep slidstep],'Value',OPTIONS.viewt0);
    set(handles.sliderSignal,'max',size(DATA.ORG,2) - (OPTIONS.viewt1 - OPTIONS.viewt0))
    set(handles.checkboxZeromean,'value',OPTIONS.zeromean);
    set(handles.zmchannels,'String','');
    
    doBaseline(handles,1);
    
end

%% --- Outputs from this function are returned to the command line.
function varargout = selectBSMbeat_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if nargout>0
    varargout{1} = handles.output;
end


%%========================================================================
function doBaseline(handles,doauto)

global DATA mlas D
global OPTIONS

minqrrs=400; % 900 for LAMAP04


DATA.BSM = DATA.ORG;
% sig1000=interp1(1:size(D.signals,1),D.signals,1:1/D.Ts:size(D.signals,1),'spline')';
% DATA.BSM = sig1000(1:size(DATA.ORG,1),1:size(DATA.ORG,2));
DATA.t = (1:size(DATA.ORG,2)) * DATA.sampT;
DATA.BSM(DATA.remove==1,:) = 0;
% lowpass at 200 Hz moving average
DATA.BSM_low200 = DATA.BSM; % DIT IS DUS ONGEFILTERD!!!

% bandpass at 0.1 to 40 Hz moving average
% BSM = lowpassma(DATA.BSM,round(1/(40*DATA.sampT))));
BSM = lowpassma(DATA.BSM,round(1/(50*DATA.sampT))); % 50Hz LP
DATA.BSMfilt = BSM - ( lowpassma(BSM, round(1/(DATA.sampT*0.1)) ) ) ;

% do zeromean
if OPTIONS.zeromean==1
    if isempty(OPTIONS.zmchannels)
        DATA.BSMfilt = zeromean(DATA.BSMfilt);
        DATA.BSM_low200 = zeromean(DATA.BSM_low200);
    else
        if min(OPTIONS.zmchannels<1) || max(OPTIONS.zmchannels)>size(DATA.BSMfilt,1)
            warning('Zeromean channels out of range. Using all valid channels');
            DATA.BSMfilt = zeromean(DATA.BSMfilt);
            DATA.BSM_low200 = zeromean(DATA.BSM_low200);
        else
            [nlds ntims]=size(DATA.BSMfilt);
            % shift to zero mean
            DATA.BSMfilt=DATA.BSMfilt-ones(nlds,1)*mean(DATA.BSMfilt(OPTIONS.zmchannels,:),1);
            [nlds ntims]=size(DATA.BSM_low200);
            % shift to zero mean
            DATA.BSM_low200=DATA.BSM_low200-ones(nlds,1)*mean(DATA.BSM_low200(OPTIONS.zmchannels,:),1);
        end
    end
end


%% find baseline points --------------------------------------------------

BSMfiltextra = BSM - ( lowpassma(BSM, round(1/(DATA.sampT*5)) ) ); % for baseline drift

rrms = rms(DATA.BSMfilt);
rrms = diffrows(rrms);

if doauto
%     drawnow
%     figure(21);
%     clf
%         plot(rrms);
%         drawnow
%         figure(1);
        fprintf('rrms mean%0.4f, std: 0.4f\n',mean(rrms),std(rrms));
    %     qrss=find(rrms>0.5*max(rrms));
    qrss=find(rrms>mean(rrms)+2*std(rrms));
    qrss(qrss < 1 )=[]; % was 300
    
    DATA.BEATS = qrss(1) - DATA.PQint;
    ignoredint=0; % ignored intervals in current rr
    for i=2:length(qrss)
        if qrss(i) - qrss(i-1)+ignoredint > minqrrs % minimal interval between qrss, was 450
            DATA.BEATS = [DATA.BEATS; qrss(i) - DATA.PQint];
            ignoredint=0;
        else
            ignoredint=ignoredint+qrss(i)-qrss(i-1); % add ignored intervals to next
        end
    end
    
    %hack for baseline drift
    if (length(DATA.BEATS) < size(DATA.BSM,2)/2000 || max(diff(DATA.BEATS))>6000) && size(BSMfiltextra,2) > 2000
        rrms = rms(BSMfiltextra);
        qrss=find(rrms>0.5*max(rrms(30:end)));
        qrss(qrss < 1 )=[]; % was 400
        DATA.BEATS = qrss(1) - DATA.PQint;
        ignoredint=0;
        for i=2:length(qrss)
            if qrss(i) - qrss(i-1) +ignoredint > minqrrs % minimal interval between qrss
                DATA.BEATS = [DATA.BEATS; qrss(i) - DATA.PQint];
                ignoredint=0;
            else
                ignoredint=ignoredint+qrss(i)-qrss(i-1);
            end
        end
    end
end
DATA.BEATS(DATA.BEATS<=0)=[];

if ~isempty(DATA.BEATS)
    zerotims=zeros(length(DATA.BEATS)+2,1);
    zerotims(1) = 1;
    zerotims(2:end-1) = DATA.BEATS;
    zerotims(end) = length(rrms);
    BASEL=spline(zerotims,DATA.BSM_low200(:,zerotims),1:size(DATA.BSM_low200,2));
    DATA.BSM_low200 = DATA.BSM_low200 - BASEL;
    BASEL=spline(zerotims,DATA.BSMfilt(:,zerotims),1:size(DATA.BSMfilt,2));
    DATA.BSMfilt    = DATA.BSMfilt - BASEL;
    
    DATA.BSMOUT=DATA.BSM; % No longer zeromean
%     DATA.BSMOUT = zeromean(DATA.BSM); % OUTPUT BY DEFINITION ZEROMEAN
    BASEL=spline(zerotims,DATA.BSMOUT (:,zerotims),1:size(DATA.BSMOUT ,2));
    DATA.BSMOUT = DATA.BSMOUT  - BASEL;
    
    DATA.BSMOUT(DATA.remove==1)=0; % oostep1: just to be sure
%     DATA.BSMOUT=zeromean(DATA.BSMOUT);
    
    
    if DATA.sampT~=1e-3
        error(' Only 1000Hz sampling implemented'); % oostep1: TBD: adapt 20ms window for other Fs
    end
    % BSMOUT50Hz: baseline from 20ms sample
    DATA.BSMOUT50Hz = DATA.BSM; % no zeromean in advance, single channel baseline drift easier to fit (50Hz common mode handled by averaging).  
    for i=1:length(DATA.BEATS)
        zerosamp(:,i)=mean(DATA.BSMOUT50Hz(:,DATA.BEATS(i)-10:DATA.BEATS(i)+9),2); % 20 ms average 
    end
    BASEL=spline(DATA.BEATS,zerosamp,1:size(DATA.BSMfilt,2)); % do not assume first and last sample 0
    DATA.BSMOUT50Hz=DATA.BSMOUT50Hz-BASEL;
    DATA.BSMOUT50Hz(DATA.remove==1)=0; % oostep1: just to be sure
%     DATA.BSMOUT50Hz=zeromean(DATA.BSMOUT50Hz);
    
    % BSMOUTSUB50Hz: subtract per-beat 20 noise template, baseline from 20ms sample
    DATA.BSMOUTSUB50Hz = DATA.BSM; % no zeromean in advance, single channel baseline drift easier to fit (50Hz common mode handled by averaging).  
    
%     DATA.BSMOUTSUB50Hz=zeromean(DATA.BSMOUTSUB50Hz); % zeromean up front if signal contains large, non-50Hz common mode distortion

    % Subtract 20ms noise template just after beat marker
    tempbeats=DATA.BEATS;
    tempbeats(end+1)=size(DATA.BSMOUTSUB50Hz,2)+1; % first sample of next virtual beat
    for i=1:length(tempbeats)-1
      DATA.BSMOUTSUB50Hz(:,tempbeats(i):tempbeats(i+1)-1)=subtracthum(DATA.BSMOUTSUB50Hz(:,tempbeats(i):tempbeats(i+1)-1));  
    end
    
    for i=1:length(DATA.BEATS)
        zerosamp(:,i)=mean(DATA.BSMOUTSUB50Hz(:,DATA.BEATS(i):DATA.BEATS(i)+19),2); % 20 ms average 
    end
    BASEL=spline(DATA.BEATS,zerosamp,1:size(DATA.BSMfilt,2)); % do not assume first and last sample 0
    DATA.BSMOUTSUB50Hz=DATA.BSMOUTSUB50Hz-BASEL;
    DATA.BSMOUTSUB50Hz(DATA.remove==1)=0; % oostep1: just to be sure
    
%     DATA.BSMOUTSUB50Hz=zeromean(DATA.BSMOUTSUB50Hz);
    
    
%     cf=get(0,'currentfigure');
%     fh=figure(345);
%     set(fh,'Name','bsmout(blue), bsmout50Hz(red) and bsmoutsub50hz(black)');
%     clf
%     ibeat=find(DATA.SELBEATS,1,'first');
%     range=DATA.BEATS(ibeat):DATA.BEATS(ibeat+1);
%     sigplot(DATA.BSMOUT(:,range),'',mlas,3,'blue'); 
%     hold on
%     sigplot(DATA.BSMOUT50Hz(:,range),'',mlas,3,'red'); 
%     sigplot(DATA.BSMOUTSUB50Hz(:,range),'',mlas,3,'k'); 
%     set(0,'currentfigure',cf);
    
else
%     DATA.BSMOUT = zeromean(DATA.BSM); % OUTPUT BY DEFINITION ZEROMEAN
    DATA.BSMOUT = DATA.BSM; % 20131024 oostep1, no longer zeromean
end
DATA.RMS     = rms(DATA.BSM_low200);
DATA.RMSfilt = rms(DATA.BSMfilt);


rrbeats=diff(DATA.BEATS);
% sprintf('RR [ms](rate [bpm]) median:%d (%d). IQR: %d-%d (%d-%d)\n',
try
    disp('RR [ms]:');
    quantile(rrbeats,[.25,.50,.75])
    disp('Rate [bpm]:');
    quantile(60000/rrbeats,[.25,.50,.75])
catch
    disp('No RR intervals measured');
end


if ~isfield(DATA,'SELBEATS') || ...
        length(DATA.SELBEATS) ~= length(DATA.BEATS)
    % resize SELBEATS, i.e. structure is cleared
    DATA.SELBEATS = zeros(size(DATA.BEATS));
    %     DATA.SELBEATS(5) = 1 % hack for ajm059
end
plotBSM(handles);


%%========================================================================
function plotBSM(handles)

global OPTIONS
global DATA

% do highpass
if get(handles.checkboxHighPass,'Value')==1
    filtertype=get(handles.popupfiltertype,'Value');
    if filtertype==1
        bsm = DATA.BSMfilt;
        rrms = DATA.RMSfilt;
    elseif filtertype==2
        if ~isempty(DATA.BEATS)
            bsm = DATA.BSMOUTSUB50Hz;
            if OPTIONS.zeromean==1
                if isempty(OPTIONS.zmchannels)
                    bsm = zeromean(bsm);
                else
                    if min(OPTIONS.zmchannels<1) || max(OPTIONS.zmchannels)>size(DATA.BSMfilt,1)
                        warning('Zeromean channels out of range. Using all valid channels');
                        bsm = zeromean(bsm);
                    else
                        [nlds ntims]=size(bsm);
                        % shift to zero mean
                        bsm=bsm-ones(nlds,1)*mean(bsm(OPTIONS.zmchannels,:),1);
                    end
                end
            end
            
            rrms = rms(bsm);
        else
            bsm = DATA.BSM_low200;
            rrms = DATA.RMS;
        end
    else
        error('illegal filtert type');
    end
else
    bsm = DATA.BSM_low200;
    rrms = DATA.RMS;
end

set(handles.axesBeats,'Visible','on');
set(handles.selectBSMFig,'CurrentAxes',handles.axesBeats);
delete(get(handles.axesBeats,'Children'));

rtop0=round(300/(1000*DATA.sampT));%find(rms(beat0)==max(rms(beat0)));
nr=0;
showsel=get(handles.ShowSelected,'Value');
for k=1:length(DATA.BEATS)-1
    if (DATA.BEATS(k) >= OPTIONS.viewt0 && DATA.BEATS(k+1) <= OPTIONS.viewt1) ||...
            (showsel==1 && DATA.SELBEATS(k)==1)
        nr=nr+1;
        beat = bsm(:,DATA.BEATS(k):DATA.BEATS(k+1));
        rtop = find(rms(beat)==max(rms(beat)));
%         if rtop<rtop0
%             beat=[zeros(size(beat,1),abs(rtop-rtop0)+1) beat]
%         elseif rtop>rtop0
%             beat=beat(:,abs(rtop-rtop0):end);
%         end
        if showsel==1 && DATA.SELBEATS(k)==1
            plot(rms(beat),'b','parent',handles.axesBeats);
        else
            plot(rms(beat),'r','parent',handles.axesBeats);
        end
        hold on
    end
    axis tight
end
% suppress stim artefact
yl=ylim;
axis manual
yl(2)=min(yl(2),1); 
ylim(yl);

axis off

%% bsm signals
set(handles.selectBSMFig,'CurrentAxes',handles.plotAxes);
delete(get(handles.plotAxes,'Children'));

nplts=size(DATA.LAY,1) - 1; % note: nplts may be smaller than nlds
raster=DATA.LAY(1,:);
axis([0 raster(1) 0 2*raster(2)] )
axis manual
% axis('off')
box off
funscal=OPTIONS.funscal;
hold on

nrastercols=DATA.LAY(1,1);

bsm  = bsm(:,OPTIONS.viewt0:OPTIONS.viewt1);
t = 0:size(bsm,2)-1;
t = t*.9/max(t) * nrastercols/(nrastercols);

for i=2:nplts+1
    j=DATA.LAY(i,3);
    yshift(j) = 2*raster(2)+1 - 2*DATA.LAY(i,2);
    xshift(j) = DATA.LAY(i,1)-1;
    if DATA.remove(j)
        plot(t(ceil(length(t)/2))+xshift(j),yshift(j),'xr','Markersize',40,'hittest','off');
    else
        plot(t+xshift(j),funscal*bsm(j,:)+yshift(j),OPTIONS.col,'linewidth',OPTIONS.linew,'hittest','off')
        plot([t(1)+xshift(j) t(end)+xshift(j)],[yshift(j) yshift(j)],':k','hittest','off')
    end
    text(xshift(j)+.4,yshift(j)+.25,num2str(j));
end

extr = extremes(bsm);
plot(t(extr(3))+xshift(extr(2)),OPTIONS.funscal*extr(1)+yshift(extr(2)),'*b');
plot(t(extr(6))+xshift(extr(5)),OPTIONS.funscal*extr(4)+yshift(extr(5)),'*r');

xl=[0.1  0.1];
yl=[0.1  1.1];
set(line(xl,yl),'color','k')
text(0.2, 0.6,sprintf('%0.3f mV',1/OPTIONS.funscal),'Tag','mVscale');
text(-DATA.LAY(1,1)/10, -.25,sprintf('%d %s %0.3f %s',size(bsm,2),' samples; ',(size(bsm,2) - 1)*DATA.sampT,'s'),'Tag','timetext');

%% RMS
set(handles.selectBSMFig,'CurrentAxes',handles.RMSaxes);
delete(get(handles.RMSaxes,'Children'));

maxA=min(max(rrms),1); % oostep1 limit y scale
% maxA=max(rrms);
% maxA=min(maxA,1*std(rrms));% suppress stim artefact

minA=0; % min(rrms); % should be 0 allways
rrms = rrms(OPTIONS.viewt0:round(OPTIONS.viewt1));
rmst = DATA.t(OPTIONS.viewt0:round(OPTIONS.viewt1));

plot(rmst,rrms,'k','hittest','off','parent',handles.RMSaxes);

hold on
beats = DATA.BEATS;
beats = beats - OPTIONS.viewt0;
beats(beats<0)=[];
beats(beats > length(rrms))=[];
for i=1:length(DATA.SELBEATS)
    if DATA.SELBEATS(i)==1 && ...
            DATA.BEATS(i) > OPTIONS.viewt0 && DATA.BEATS(i+1) < OPTIONS.viewt1
        plot(DATA.t(DATA.BEATS(i):DATA.BEATS(i+1)),zeros(size(DATA.t(DATA.BEATS(i):DATA.BEATS(i+1)))),'r','linewidth',10,'parent',handles.RMSaxes);
        hold on;
    end
end

plot(rmst(beats),rrms(beats),'or','Markersize',10,'parent',handles.RMSaxes);
axis([rmst(1) rmst(end) minA maxA]);

% Overlay RMS with signal
os=get(handles.OverlaySignal,'Value');
siglist=get(handles.OverlaySignal,'String');
if os>1 % 1 is none
    sig=os-1;;
    if sig<=size(DATA.ORG,1)
        osig=DATA.ORG(sig,:);
    else
        osig=DATA.EGM(sig-size(DATA.ORG,1),:);
    end
    osig=osig- ( lowpassma(osig,round(1/(DATA.sampT *0.1)) ) );
    osig=osig-mean(osig);
    osig=(minA+maxA)/2+osig/(2*std(osig))*(maxA-minA);
    osigs=osig(OPTIONS.viewt0:round(OPTIONS.viewt1));
    plot(rmst,osigs,'r');
end

set(handles.nrselected,'String',num2str(sum(DATA.SELBEATS)));

set(handles.RMSaxes ,'ButtonDownFcn','selectBSMbeat(''RMSaxes_ButtonDownFcn'' ,gcbo,[],guidata(gcbo))');
set(handles.plotAxes,'ButtonDownFcn','selectBSMbeat(''plotAxes_ButtonDownFcn'' ,gcbo,[],guidata(gcbo))');


%% --- Executes on slider movement.
function sliderSignal_Callback(hObject, eventdata, handles)
% hObject    handle to sliderSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global DATA
global OPTIONS

OPTIONS.viewt0 = round(get(handles.sliderSignal,'Value'));
steps = get(handles.sliderSignal,'sliderstep');
OPTIONS.viewt1 = min(size(DATA.ORG,2),round(OPTIONS.viewt0 + (steps(2) * size(DATA.ORG,2))));
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

global OPTIONS
global DATA
OPTIONS.viewt0 = round(get(handles.sliderSignal,'Value')); % PO was commented out
steps = get(handles.sliderSignal,'sliderstep');
steps = steps / 2;
set(handles.sliderSignal,'sliderstep',steps);

% OPTIONS.viewt0 = round(get(handles.sliderSignal,'Value') + (steps(2)/2 * size(DATA.ORG,2))); % PO
OPTIONS.viewt1 = min(size(DATA.ORG,2),round(OPTIONS.viewt0 + (steps(2) * size(DATA.ORG,2))));
set(handles.sliderSignal,'max',size(DATA.ORG,2) - (OPTIONS.viewt1 - OPTIONS.viewt0))
plotBSM(handles)

%% --- Executes on button press in pushZoomOut.
function pushZoomOut_Callback(hObject, eventdata, handles)
% hObject    handle to pushZoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global DATA
global OPTIONS

OPTIONS.viewt0 = round(get(handles.sliderSignal,'Value'));  % PO was commented out
steps = get(handles.sliderSignal,'sliderstep');
steps = steps * 2;
set(handles.sliderSignal,'sliderstep',steps);
% OPTIONS.viewt0 =max(0,round(get(handles.sliderSignal,'Value')-steps(2)/2*size(DATA.ORG,2)));%PO
OPTIONS.viewt1 = min(size(DATA.ORG,2),round(OPTIONS.viewt0 + (steps(2) * size(DATA.ORG,2))));
set(handles.sliderSignal,'max',size(DATA.ORG,2) - (OPTIONS.viewt1 - OPTIONS.viewt0))

plotBSM(handles)

%% --- Executes on mouse press over axes background.
function RMSaxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to RMSaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DATA
global OPTIONS
point1=get(hObject,'currentpoint');
selTime=round(point1(1)/DATA.sampT);
selMarker = find(abs(DATA.BEATS-selTime) < 20);
maxA=max(DATA.RMS);
if strcmp(get(handles.selectBSMFig,'SelectionType'),'extend')% add zeropoint
    if isempty(find(abs(DATA.BEATS-selTime) < 300)) % no marker found within 300 ms
        a = find(find(DATA.BEATS-selTime < 0));
        if isempty(a)
            DATA.BEATS = [selTime;DATA.BEATS];
            DATA.SELBEATS = [0 ; DATA.SELBEATS];
        else
            a=a(end);
            if a ~= length(DATA.BEATS)
                DATA.BEATS = [DATA.BEATS(1:a); selTime; DATA.BEATS(a+1:end)];
                DATA.SELBEATS = [DATA.SELBEATS(1:a); 0 ; DATA.SELBEATS(a+1:end)];
            else
                DATA.BEATS = [DATA.BEATS(1:a); selTime; ];
                DATA.SELBEATS = [DATA.SELBEATS(1:a); 0  ];
            end
        end
    end
elseif strcmp(get(handles.selectBSMFig,'SelectionType'),'alt')% select point
    if size(DATA.SELBEATS,1) ~= size(DATA.BEATS,1)
        DATA.SELBEATS = zeros(size(DATA.BEATS));
    end
    markers = DATA.BEATS-selTime;
    markers(markers<0) =10000;
    selBeat = find(markers == min(markers))-1;
    if selBeat > 1
        DATA.SELBEATS(selBeat) = ~DATA.SELBEATS(selBeat);
    end
elseif ~isempty(selMarker)
    figpoint1=get(handles.selectBSMFig,'currentpoint');
    rect=[figpoint1(1,1), figpoint1(1,2), 2, maxA];
    rbbox;
    point1=get(hObject,'currentpoint');
    DATA.BEATS(selMarker) = round(point1(1)/DATA.sampT);
    doBaseline(handles,0);
else
    
    %     delete(findobj('Tag','pline1')); delete(findobj('Tag','pline2'));
    %     line([point1(1,1) point1(1,1)],[0 maxA],'Color','r','parent',hObject,'Tag','pline1','hittest','off');
    %     rbbox;
    %
    %     % return figure units
    %     point2 = get(hObject,'CurrentPoint');% button up detected
    %
    %     line([point2(1,1) point2(1,1)],[0 maxA],'Color','r','parent',hObject,'Tag','pline2','hittest','off');
    %     xmax=max(point1(1,1),point2(1,1));
    %     xmin=min(point1(1,1),point2(1,1));
    %     if xmax < size(DATA.BSM,2) - 1 %add one second min width
    %         xmax = max(xmin+1,xmax);
    %     else
    %         xmin = min(xmin,xmax-1);
    %     end
    % %     disp(num2str([xmin xmax],4));
    %     OPTIONS.viewt0 = max(1,round(xmin/DATA.sampT));
    %     OPTIONS.viewt1 = min(size(DATA.BSM,2),round(xmax/DATA.sampT));
    %     delete(findobj('Tag','pline1')); delete(findobj('Tag','pline2'));
    %     updateSignalSlider(handles);
    %     set(handles.sliderSignal,'Value',(OPTIONS.viewt0+OPTIONS.viewt1)/2);
end

plotBSM(handles)



% --- Executes on mouse press over axes background.
function plotAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to plotAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DATA

% 	yshift(j)=2*raster(2)+1-2*DATA.LAY(i,2);
% 	xshift(j)=DATA.LAY(i,1)-1;

point1=get(hObject,'currentpoint');
x=floor(1+point1(1,1));
y=-round((point1(1,2)-DATA.LAY(1,2)*2-1)/2);
% disp(num2str([point1(1,1) point1(1,2) x y]))
a=find(DATA.LAY(:,1)==x & DATA.LAY(:,2)==y);
if ~isempty(a)
    DATA.remove(DATA.LAY(a,3))=xor(DATA.remove(DATA.LAY(a,3)),1);
end
doBaseline(handles,0);

% --- Executes on button press in checkboxHighPass.
function checkboxHighPass_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxHighPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxHighPass

% global DATA

% DATA.BSM=DATA.ORG;
% if get(hObject,'Value')
% 	DATA.BSM=DATA.BSM-lowpassma(DATA.BSM,round(DATA.len20ms/DATA.sampT)); %?50Hz
% end
plotBSM(handles);

% --- Executes during object creation, after setting all properties.
function sliderSignal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in checkboxZeromean.
function checkboxZeromean_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxZeromean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global OPTIONS

OPTIONS.zeromean=get(hObject,'value');
doBaseline(handles,0);


% --- Executes on scroll wheel click while the figure is in focus.
function selectBSMFig_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to selectBSMFig (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)

global OPTIONS

OPTIONS.funscal = max(0.3,OPTIONS.funscal - eventdata.VerticalScrollCount/5.0);
plotBSM(handles);

% --- Executes when user attempts to close selectBSMFig.
function selectBSMFig_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to selectBSMFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);

% --- Executes during object creation, after setting all properties.
function selectBSMFig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectBSMFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes during object creation, after setting all properties.
function recompBaseline_Button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recompBaseline_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uipanelControls_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanelControls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in recompBaseline_Button.
function recompBaseline_Button_Callback(hObject, eventdata, handles)
% hObject    handle to recompBaseline_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

doBaseline(handles,1);


function PQint_Callback(hObject, eventdata, handles)
% hObject    handle to PQint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PQint as text
%        str2double(get(hObject,'String')) returns contents of PQint as a double
global DATA
pqintms=str2double(get(hObject,'String'));
pqintsamp = round(pqintms/1000/DATA.sampT);
DATA.BEATS=DATA.BEATS+DATA.PQint-pqintsamp;
DATA.PQint=pqintsamp;
doBaseline(handles,0);








% --- Executes during object creation, after setting all properties.
function PQint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PQint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zmchannels_Callback(hObject, eventdata, handles)
% hObject    handle to zmchannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zmchannels as text
%        str2double(get(hObject,'String')) returns contents of zmchannels as a double
global OPTIONS
chstr=get(hObject,'String');
chstr=regexprep(chstr,'[^0-9:\s,]',''); % remove illegal characters
set(hObject,'String',chstr);
try
    OPTIONS.zmchannels=eval(['[' chstr ']']);
catch
    OPTIONS.zmchannels=[];
    set(hObject,'String','');
end
if OPTIONS.zeromean
    doBaseline(handles,0);
end


% --- Executes during object creation, after setting all properties.
function zmchannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zmchannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showonthorax.
function showonthorax_Callback(hObject, eventdata, handles)
% hObject    handle to showonthorax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global OPTIONS
global DATA
% do highpass
if get(handles.checkboxHighPass,'Value')==1
        filtertype=get(handles.popupfiltertype,'Value');
    if filtertype==1
        bsm = DATA.BSMfilt;
    else % For now only one oter option: 2=subtracthum
        
        if OPTIONS.zeromean==1
            if isempty(OPTIONS.zmchannels)
                bsm = zeromean(DATA.BSMOUTSUB50Hz);
            else
                if min(OPTIONS.zmchannels<1) || max(OPTIONS.zmchannels)>size(DATA.BSMfilt,1)
                    warning('Zeromean channels out of range. Using all valid channels');
                    bsm = zeromean(DATA.BSMOUTSUB50Hz);
                else
                    [nlds ntims]=size(DATA.BSMOUTSUB50Hz);
                    % shift to zero mean
                    bsm=DATA.BSMOUTSUB50Hz-ones(nlds,1)*mean(DATA.BSMOUTSUB50Hz(OPTIONS.zmchannels,:),1);
                end
            end
        end
        
    end
    
else
    bsm = DATA.BSM_low200;
end
bsm=lowpassma(bsm,3);
ShowBSMonThorax(bsm(:,OPTIONS.viewt0:OPTIONS.viewt1),OPTIONS.leadsys,find(DATA.remove==0)); % plot zeromean data on thorax
% ShowBSMonThorax(DATA.RMS(OPTIONS.viewt0:OPTIONS.viewt1),OPTIONS.leadsys); % plot zeromean data on thorax



% --- Executes on button press in readmatbutton.
function readmatbutton_Callback(hObject, eventdata, handles)
% hObject    handle to readmatbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
persistent defpath
global DATA
global OPTIONS
global D mlas
[OPTIONS.FileName,OPTIONS.PathName,FilterIndex] = uigetfile({'*.mat','*.*'},'Select data file.',defpath);
defpath=OPTIONS.PathName;
disp(OPTIONS.FileName);

load(fullfile(OPTIONS.PathName,OPTIONS.FileName),'DATA','D','mlas');

% Haal EGM uit D (resample) als geen veld EGM
egmrange=65:79;
if ~isfield(DATA,'EGM')
    if size(D.signals,2)<max(egmrange)
        comment=D.comment;
        filepath=D.filepath;
        D=Readbdf(D.filepath,1:max(egmrange),D.startpos,D.endpos,'noui');
        D.comment=comment;
        D.filepath=filepath;
    end
    DATA.EGM=interp1(1:size(D.signals,1),D.signals(:,egmrange),1:1/D.Ts:size(D.signals,1),'cubic')';
end


set(handles.PQint,'String',num2str(DATA.PQint*DATA.sampT*1000));
OPTIONS.viewt0=1;
OPTIONS.viewt1=ceil(max(100,size(DATA.ORG,2)/4));
set(handles.sliderSignal,'sliderStep',[0.01,0.1])
set(handles.sliderSignal,'min',1);


OPTIONS.zeromean=0;
OPTIONS.zmchannels=[];
set(handles.checkboxHighPass,'value',0);

set(handles.sliderSignal,'value',1);
slidstep = ((OPTIONS.viewt1 - OPTIONS.viewt0) / size(DATA.ORG,2));
set(handles.sliderSignal,'Sliderstep',[0.1*slidstep slidstep],'Value',OPTIONS.viewt0);
set(handles.sliderSignal,'max',size(DATA.ORG,2) - (OPTIONS.viewt1 - OPTIONS.viewt0))
set(handles.checkboxZeromean,'value',OPTIONS.zeromean);
set(handles.zmchannels,'String','');

%Init OverlaySignal list
siglist=[{'(none)'};cellstr(num2str((1:size(DATA.ORG,1))'))];
if isfield(DATA,'EGM');
    for i=1:size(DATA.EGM,1)
        siglist{end+1}=sprintf('egm%03d',i);
    end
end
set(handles.OverlaySignal,'String',siglist)
set(handles.OverlaySignal,'Value',1);

set(handles.Exportbutton,'Visible','on');

% plotBSM(handles);
doBaseline(handles,0);


% --- Executes on button press in ShowSelected.
function ShowSelected_Callback(hObject, eventdata, handles)
% hObject    handle to ShowSelected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ShowSelected
plotBSM(handles);


% --- Executes on selection change in OverlaySignal.
function OverlaySignal_Callback(hObject, eventdata, handles)
% hObject    handle to OverlaySignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns OverlaySignal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from OverlaySignal

plotBSM(handles);



% --- Executes during object creation, after setting all properties.
function OverlaySignal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OverlaySignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Exportbutton.
function Exportbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Exportbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global OPTIONS
global DATA D mlas
[FileName,PathName,FilterIndex] = uiputfile({'*.mat','*.*'},'Save as.',fullfile(OPTIONS.PathName,OPTIONS.FileName));
if PathName~=0
    OPTIONS.FileName=FileName;
    OPTIONS.PathName=PathName;
    save(fullfile(PathName,FileName),'DATA','D','mlas');
end



% --- Executes on selection change in popupfiltertype.
function popupfiltertype_Callback(hObject, eventdata, handles)
% hObject    handle to popupfiltertype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupfiltertype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupfiltertype
plotBSM(handles);

% --- Executes during object creation, after setting all properties.
function popupfiltertype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupfiltertype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
