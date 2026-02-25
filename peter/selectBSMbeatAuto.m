function varargout = selectBSMbeatAuto(varargin)
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

% oostep1:  DATA.PQint set to 120 ms instead of 180 samples.
%           round lowpassma window
%           Add ui control for DATA.PQint
%           PQint init at 220 ms, at 280 ms
% 20120622 oostep1: select channels for zeromean, ShowBSMonThorax, doBaseline when changing PQInt
% 20120813 oostep1: zoom to center.


% Last Modified by GUIDE v2.5 18-Dec-2012 17:58:50

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

% global handles
set(handles.selectBSMFig,'Toolbar','figure')
set(handles.selectBSMFig,'menubar','figure')

if isfield( DATA, 'BEATS') 
    DATA = rmfield(DATA,{'BEATS','PWAVES','BSMraw','BSM','t'});
end


OPTIONS.funscal=0.5;
OPTIONS.zeromean=0;
OPTIONS.col='blue';
OPTIONS.linew=1;
OPTIONS.label='';
OPTIONS.leadsys='ams';%'pigbsm'/'ams'; % default leadsystem for ShowBSMonThorax
OPTIONS.zmchannels=[];
OPTIONS.selectLead = -1;

if length(varargin) < 1
    % 	error('This routine needs at least two parameters');
    readmatbutton_Callback(handles.readmatbutton,eventdata,handles); % read DATA from file
else
    set(handles.Exportbutton,'Visible','off'); % export handled by caller
    DATA.ORG = varargin{1};
    DATA.sampT=1/1000;
    zoomfull = 0;
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
                case 'zoomfull'
                    zoomfull=varargin{pp+1};pp=pp+2;
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
    
    % BSM = lowpassma(DATA.ORG,1/(40*DATA.sampT));
    sumData = mean(abs(DATA.ORG - ( lowpassma(DATA.ORG, 2 / ( DATA.sampT ) ) )),2);
    % sumData = mean((DATA.ORG),2);
    DATA.remove = sumData > 7 * std(sumData);
    pqint=170; % PQ interval in ms
    DATA.PQint = round(pqint/1000/DATA.sampT); % default PQinterval + a part of the QRS interval, in samples
    set(handles.PQint,'String',num2str(pqint));
    
    OPTIONS.viewt0=1;
    %     OPTIONS.viewt1=ceil(max(100,size(DATA.ORG,2)/4)); % PO
    if zoomfull
        OPTIONS.viewt1=size(DATA.ORG,2);
    else
        OPTIONS.viewt1=ceil(min(5000,size(DATA.ORG,2)/2));
    end
    
      
    % Choose default command line output for selectBSMbeat
    handles.output = hObject;
    
    % Update handles structure
    guidata(hObject, handles);
    
    % UIWAIT makes selectBSMbeat wait for user response (see UIRESUME)
    % uiwait(handles.selectBSMFig);
    
    set(handles.sliderSignal,'sliderStep',[0.01,0.1])
    set(handles.sliderSignal,'min',0);  
    set(handles.sliderSignal,'value',0);
    set(handles.sliderSignal,'max',0.9); 
    
    set(handles.checkboxZeromean,'value',OPTIONS.zeromean);
    set(handles.zmchannels,'String','');
    
    doBaseline(handles,1);
    
    %Init OverlaySignal list
    siglist=[{'(none)'};cellstr(num2str((1:size(DATA.ORG,1))'))];
    set(handles.OverlaySignal,'String',siglist)
    set(handles.OverlaySignal,'Value',1);
    
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

global DATA
global OPTIONS

DATA.BSM = DATA.ORG;
DATA.t = (1:size(DATA.ORG,2)) * DATA.sampT;
DATA.BSM(DATA.remove==1,:) = 0;
DATA.BSMraw = zeromean(DATA.BSM); 


% do zeromean
if ( ( isfield(OPTIONS,'prevzeromean') && OPTIONS.prevzeromean == OPTIONS.zeromean && OPTIONS.zeromean == 1 ) || ...
       OPTIONS.zeromean == 1 )
    disp('do zeromean');
    % bandpass at 0.1 to 40 Hz moving average
    BSM = lowpassma(DATA.BSM,round( 1 / (50*DATA.sampT))); % 50Hz LP
    DATA.BSMfilt = BSM - ( lowpassma(BSM, 1 / ( 0.1 * DATA.sampT ) ) );

    if isempty(OPTIONS.zmchannels)
        DATA.BSMfilt = zeromean(DATA.BSMfilt);
    else
        if min(OPTIONS.zmchannels<1) || max(OPTIONS.zmchannels)>size(DATA.BSMfilt,1)
            warning('Zeromean channels out of range. Used all valid channels');
            DATA.BSMfilt = (DATA.BSMfilt);
        else
            [nlds, ~]=size(DATA.BSMfilt);
            % shift to zero mean
            DATA.BSMfilt=DATA.BSMfilt-ones(nlds,1)*mean(DATA.BSMfilt(OPTIONS.zmchannels,:),1);
        end
    end
elseif ~isfield(DATA,'BSMfilt') || size(DATA.BSMfilt,2) ~= size(DATA.ORG,2)
    % bandpass at 0.1 to 40 Hz moving average
    BSM = lowpassma(DATA.BSM,round( 1 / (50*DATA.sampT))); % 50Hz LP
    DATA.BSMfilt = BSM - ( lowpassma(BSM, 1 / ( 0.1 * DATA.sampT ) ) );
    disp('compute filtered signals');
end

OPTIONS.prevzeromean = OPTIONS.zeromean;

%% find baseline points --------------------------------------------------

% BSMfiltextra = DATA.BSM  - ( lowpassma(DATA.BSM , round(1 / ( 10.0 * DATA.sampT ) ) ) ); % for baseline drift


if doauto
    disp('determine beats');
    qrrms = rms(DATA.BSMfilt);
       
%     qrss = find(qrrms > mean(qrrms) + 5 * std(qrrms) );
    maxvalue=max(qrrms(1:1000));%mean(qrrms) + 4 * std(qrrms);
    maxValueOrg = maxvalue;
    qrsPeaki=find(qrrms(1:1000)==maxvalue);
    qrss=[];
    i=qrsPeaki-1;
    maxvalue= qrrms(i) * 0.9;
    while i < length(qrrms)
        maxvalue = maxvalue -0.0001;%mean(qrrms) + 1 * std(qrrms);
        if abs(maxvalue - qrrms(i) ) > 5 * maxValueOrg
            maxvalue = maxValueOrg*0.8;
        end
        if qrrms(i) > maxvalue
            maxvalue = qrrms(i);
            qrsPeaki = i;
        elseif i > qrsPeaki + 100
            qrss=[qrss qrsPeaki-DATA.PQint];
            i = i + DATA.PQint;
            qrsPeaki = length(qrrms);
%             maxvalue = mean(qrrms) + 1 * std(qrrms);
        end
         
        i= i+1;
    end
    
    if qrss(1) < 300
        qrss(1) =[];
    end
    doShow =0;
    rrms = rms(lowpassma(DATA.BSM,21));
%     rrms = rms(DATA.BSM);
    DATA.BEATS = [];%qrss(1);
    DATA.PWAVES=[];
    DATA.PEAKTWAVES=[];
    DATA.JPOINTS=[];
    Ponset =-1;
    
    
    for i=1:length(qrss) -1
        
        
%         if qrss(i) - qrss(i-1) > 190 % minimal interval between qrss
            peakQrs = min(190,qrss(i)-1);%,floor((qrss(i) - qrss(i-1))/3));       
            
            if i < length(qrss)
                sig = rrms(max(1,qrss(i)-peakQrs):min(length(rrms),find(rrms == min(rrms(qrss(i+1)-peakQrs:qrss(i+1))))));    
            else
                sig = rrms(qrss(i)-peakQrs:min(length(rrms),qrss(i) + 250));
                sig = sig(:,find(sig == min(sig(1:peakQrs-80))):find(sig == min(sig(peakQrs+1:end))));
                
            end
            b=find(sig==min(sig(length(sig)-peakQrs:end)));
            sig = baselinecor(sig,find(sig(1:peakQrs)== min(sig(1:peakQrs))),b(end));
            dsig= abs(diffrows(sig));
            
           
            delta = 3;

            minval=sig(peakQrs);
            dminval=dsig(peakQrs);
            iminval = find(dsig(1:peakQrs)==max(dsig(max(1,peakQrs-30):peakQrs)));
            prevNoiseVal =1;
            keepiMinVal=1;
            for j = iminval:-1:max(50,iminval-100)                
%                 sigtest = baselinecor(sig,j,length(sig),1); 
%                 meanVal = mean(abs(sigtest(j-delta:j+delta) ));
                dmeanVal = mean(abs(dsig(j-delta:j+delta)  ));
                meanVal = mean(abs(sig(j-delta:j+delta)  ));
                sigtest = baselinecor(sig,j,length(sig),1); 
                noiseVal = mean(abs(sigtest(j-20:j)));
                
                if noiseVal  > mean(abs(sigtest(j:j+4))) % (dmeanVal < dminval && iminval - j < 25)||   (meanVal < minval && iminval - j < 10 )
                    if doShow
                        figure(1);clf;plot(sig);hold on; plot(sigtest,'r');plot(j,sigtest(j),'ko');plot(dsig*10,'g');drawnow
                        disp(num2str([peakQrs-j meanVal noiseVal ]))
                    end
                    minval = min(minval,meanVal);
                    dminval = min(dminval,dmeanVal);
                    iminval = j;     
                else
                    stop=1;
                    if noiseVal < prevNoiseVal
                        keepiMinVal = iminval;
                    end                    
                end
                prevNoiseVal = min(noiseVal,prevNoiseVal);

            end
            
            iminval = keepiMinVal;

            onset = peakQrs - iminval;
            
            % Pwave detection            
            maxval = sig(1:peakQrs - onset);
            imaxval = peakQrs - onset;
            founddecent=0;
            for j = peakQrs - onset-50:-1: max(1,peakQrs - onset - 90) %peakQrs - onset - 30
                if sig(j) < sig(j+1) && sig(j) < sig(j+2) && sig(j) < sig(j+3) && sig(j) < sig(j+4)
                    maxval = sig(j);
                    founddecent =1;
                    imaxval = j;        
                elseif founddecent
                    break;
                end
            end 
            imaxval =  find(sig== max(sig(1:imaxval(1))));
            iminval = imaxval;
            minval = 1000;%sig(iminval);
            for j = imaxval:-1: max(delta+1,imaxval - 200) %peakQrs - onset - 30
                meanVal = mean(abs(sig(j-delta:j+delta)  ));
                if meanVal < minval  && iminval - j < 10 
                    if doShow
                        sigtest = baselinecor(sig,j,length(sig),1); 
                        figure(2);clf;plot(sig);hold on; plot(sigtest,'r');plot(j,sigtest(j),'ko');drawnow
%                         disp(num2str([peakQrs-j meanVal  ]))
                    end

                    minval = meanVal;
                    iminval = j;                                      
                end
            end
            Ponset = peakQrs - iminval;

            if iminval > delta + 10
%                 Ponset = onset;
                DATA.PWAVES=[DATA.PWAVES ; qrss(i) - Ponset ];               
                pwave =1;
            else
                Ponset = onset;
                pwave =0;
            end
            if ~pwave  &&0               
                iminval = peakQrs-11;
                minval=sig(iminval);
                for j = iminval:-1:10                
                    meanVal = mean(abs(sig(j-delta:j+delta)  ));
                    if meanVal < minval  && iminval - j < 10 
                        if doShow
                            sigtest = baselinecor(sig,j,length(sig),1); 
                            figure(1);clf;plot(sig);hold on; plot(sigtest,'r');plot(j,sigtest(j),'ko');drawnow
                            disp(num2str([peakQrs-j meanVal  ]))
                        end
                        minval = meanVal;
                        iminval = j;                                      
                    end

                end
                onset = peakQrs - iminval; 
            end
          
            
            if Ponset ~= onset && onset > 170 && abs(onset -Ponset) < 100
            
                minval=sig(peakQrs);
                iminval = peakQrs-11;
                for j = iminval:-1:Ponset-100
                
%                 meanVal = mean(abs(sigtest(j-delta:j+delta) ));
                    meanVal = mean(abs(sig(j-delta:j+delta)  ));
                    if meanVal < minval  && iminval - j < 10 
                        if doShow
                            sigtest = baselinecor(sig,j,length(sig),1); 
                            figure(1);clf;plot(sig);hold on; plot(sigtest,'r');plot(j,sigtest(j),'ko');drawnow
%                             disp(num2str([peakQrs-j meanVal  ]))
                        end


                        minval = meanVal;
                        iminval = j;                                      
                    end
                end
                disp(num2str([onset peakQrs - iminval]))
                onset = peakQrs - iminval;
            
            end
            
            
            
%             DATA.BEATS = [DATA.BEATS; qrss(i) - Ponset ];
            DATA.BEATS = [DATA.BEATS; qrss(i) - onset ];
            if length(DATA.BEATS) > 1
%                 sig = lowpassma(rrms(DATA.BEATS(end-1):qrss(i) - Ponset ),21);
               
                
                iminval = 80;
                minval=sig(iminval);                
                for j = iminval : min(250 ,length(sig)-delta)
                
%                 meanVal = mean(abs(sigtest(j-delta:j+delta) ));
                    meanVal = mean(abs(sig(j-delta:j+delta)  ));
                    if meanVal < minval  && iminval - j < 10 
                        if doShow
                            sigtest = baselinecor(sig,j,length(sig),1); 
                            figure(3);clf;plot(sig);hold on; plot(sigtest,'r');plot(j,sigtest(j),'ko');drawnow
%                             disp(num2str([peakQrs-j meanVal  ]))
                        end


                        minval = meanVal;
                        iminval = j;                                      
                    end

                end
                
                DATA.PEAKTWAVES =   [DATA.PEAKTWAVES; DATA.BEATS(end-1) + find(sig == max(sig(iminval:end))) ];
                
                DATA.JPOINTS = [DATA.JPOINTS; DATA.BEATS(end-1) + iminval ];
            end

%         end
    end      
end
% DATA.BEATS = DATA.PWAVES
DATA.BEATS(DATA.BEATS<=0)=[];
tic
if ~isempty(DATA.BEATS)
    if DATA.BEATS(1)==1
        zerotims=zeros(length(DATA.BEATS)+1,1);
        zerotims(1:end-1) = DATA.BEATS;
        zerotims(end) = length(rrms);
    else
        zerotims=zeros(length(DATA.BEATS)+2,1);
        zerotims(1) = 1;
        zerotims(2:end-1) = DATA.BEATS;
        zerotims(end) = length(rrms);       
    end
    disp('baseline')
    tic 
%     BASEL=spline(zerotims,DATA.BSMraw(:,zerotims),1:size(DATA.BSMraw,2));
    DATA.BSMraw = DATA.BSMraw - spline(zerotims,DATA.BSMraw(:,zerotims),1:size(DATA.BSMraw,2));
%     BASEL=spline(zerotims,DATA.BSMfilt(:,zerotims),1:size(DATA.BSMfilt,2));
    DATA.BSMfilt = DATA.BSMfilt - spline(zerotims,DATA.BSMfilt(:,zerotims),1:size(DATA.BSMfilt,2));    
    DATA.BSMOUT  = DATA.BSMraw;%zeromean(DATA.BSMraw); % OUTPUT BY DEFINITION ZEROMEAN    
else
    DATA.BSMOUT = DATA.BSMraw;%zeromean(DATA.BSMfilt); % OUTPUT BY DEFINITION ZEROMEAN
end
DATA.RMS     = rms(DATA.BSMraw);
DATA.RMSfilt = rms(DATA.BSMfilt);
toc

% rrbeats=diff(DATA.BEATS);
% sprintf('RR [ms](rate [bpm]) median:%d (%d). IQR: %d-%d (%d-%d)\n',
% disp('RR [ms]:');
% quantile(rrbeats,[.25,.50,.75])
% disp('Rate [bpm]:');
% quantile(60000/rrbeats,[.25,.50,.75])



if ~isfield(DATA,'SELBEATS') || ...
        length(DATA.SELBEATS) ~= length(DATA.BEATS)
    % resize SELBEATS, i.e. structure is cleared
    DATA.SELBEATS = zeros(size(DATA.BEATS));
end
tic
plotBSM(handles);
toc

%%========================================================================
function plotBSM(handles)

global OPTIONS
global DATA

% OPTIONS.viewt0=5000;
% OPTIONS.viewt1=10000;

% do highpass
if get(handles.checkboxHighPass,'Value')==1
    bsm = DATA.BSMfilt;
    rrms = DATA.RMSfilt;
else
    bsm = DATA.BSMraw;
    rrms = DATA.RMS;
end

% set(handles.axesBeats,'Visible','on');
set(handles.selectBSMFig,'CurrentAxes',handles.axesBeats);
delete(get(handles.axesBeats,'Children'));

rtop0=round(300/(1000*DATA.sampT));%find(rms(beat0)==max(rms(beat0)));
nr=0;
maxrms= min(max(rrms),2);
showsel=get(handles.ShowSelected,'Value');
for k=1:length(DATA.BEATS)-1
    if (DATA.BEATS(k) >= OPTIONS.viewt0 && DATA.BEATS(k+1) <= OPTIONS.viewt1) ||...
            (showsel==1 && DATA.SELBEATS(k)==1)
        nr=nr+1;
        beat = bsm(:,DATA.BEATS(k):DATA.BEATS(k+1));
        rtop = find(rms(beat)==max(rms(beat)));
        if rtop<rtop0
            beat=[zeros(size(beat,1),abs(rtop-rtop0)+1) beat];
        elseif rtop>rtop0
            beat=beat(:,abs(rtop-rtop0):end);
        end

        if showsel==1 && DATA.SELBEATS(k)==1
            plot(rms(beat)/maxrms,'b','parent',handles.axesBeats);
        else
            plot(rms(beat)/maxrms,'r','parent',handles.axesBeats);
        end
        hold on
    end
    
end
axis([0 1200 0 1]);
% suppress stim artefact
yl=ylim;
yl(2)=min(yl(2),1);
ylim(yl);


%% bsm signals
set(handles.selectBSMFig,'CurrentAxes',handles.plotAxes);
delete(get(handles.plotAxes,'Children'));

nplts=size(DATA.LAY,1) - 1; % note: nplts may be smaller than nlds
raster=DATA.LAY(1,:);
axis([0 raster(1) 0 2 * raster(2)] )
funscal=OPTIONS.funscal;
hold on

nrastercols=DATA.LAY(1,1);

bsm  = bsm(:,OPTIONS.viewt0:OPTIONS.viewt1);
t = 0:size(bsm,2)-1;
t = t / max(t);
if raster(1) > 1
    t = t*0.9;
end
    

for i=2:nplts+1
    j=DATA.LAY(i,3);
    yshift(j) = 2 * raster(2) - 2 * (DATA.LAY(i,2) - 0.5);
    xshift(j) = DATA.LAY(i,1) - 1;
    if DATA.remove(j)
        plot(t(ceil(length(t)/2))+xshift(j),yshift(j),'xr','Markersize',40,'hittest','off');
    else
        plot(t+xshift(j),funscal*bsm(j,:)+yshift(j),OPTIONS.col,'linewidth',OPTIONS.linew,'hittest','off')
        plot([t(1)+xshift(j) t(end)+xshift(j)],[yshift(j) yshift(j)],':k','hittest','off')
    end
    text(xshift(j)+.4,yshift(j)+.25,num2str(j));
end

% extr = extremes(bsm);
% plot(t(extr(3))+xshift(extr(2)),OPTIONS.funscal*extr(1)+yshift(extr(2)),'*b');
% plot(t(extr(6))+xshift(extr(5)),OPTIONS.funscal*extr(4)+yshift(extr(5)),'*r');

% xl=[0.1  0.1];
% yl=[0.1  1.1];
% set(line(xl,yl),'color','k')
% text(0.2, 0.6,sprintf('%0.3f mV',1/OPTIONS.funscal),'Tag','mVscale');
% text(-DATA.LAY(1,1)/10, -.25,sprintf('%d %s %0.3f %s',size(bsm,2),' samples; ',(size(bsm,2) - 1)*DATA.sampT,'s'),'Tag','timetext');

%% RMS

set(handles.selectBSMFig,'CurrentAxes',handles.RMSaxes);
delete(get(handles.RMSaxes,'Children'));
rrms = rrms(OPTIONS.viewt0:round(OPTIONS.viewt1));
rmst = DATA.t(OPTIONS.viewt0:round(OPTIONS.viewt1));

beats = DATA.BEATS;
beats = beats - OPTIONS.viewt0;
beats(beats<=0)=[];
beats(beats > length(rrms))=[];

pwaves = DATA.PWAVES;
pwaves = pwaves - OPTIONS.viewt0;
pwaves(pwaves<=0)=[];
pwaves(pwaves > length(rrms))=[];
disp([ 'pwaves ' num2str(pwaves')])

twaves = DATA.PEAKTWAVES;
twaves = twaves - OPTIONS.viewt0;
twaves(twaves<=0)=[];
twaves(twaves > length(rrms))=[];

Jpnt = DATA.JPOINTS;
Jpnt = Jpnt - OPTIONS.viewt0;
Jpnt(Jpnt<=0)=[];
Jpnt(Jpnt > length(rrms))=[];

maxA = max(rrms);
% maxA = min(maxA,2);% suppress stim artefact oostp1
minA = min(rrms);
    
if OPTIONS.selectLead < 0
    plot(rmst,rrms,'k','hittest','off','parent',handles.RMSaxes);
 
    hold all
    for i=1:length(DATA.SELBEATS)
        if DATA.SELBEATS(i)==1 && ...
                DATA.BEATS(i) > OPTIONS.viewt0 && DATA.BEATS(i+1) < OPTIONS.viewt1
            plot(DATA.t(DATA.BEATS(i):DATA.BEATS(i+1)),zeros(size(DATA.t(DATA.BEATS(i):DATA.BEATS(i+1)))),'r','linewidth',10,'parent',handles.RMSaxes);
            hold on;
        end
    end
    plot(rmst(beats),rrms(beats),'or','Markersize',10,'parent',handles.RMSaxes);
    
    plot(rmst(pwaves),rrms(pwaves),'gx','Markersize',10,'parent',handles.RMSaxes);
    plot(rmst(twaves),rrms(twaves),'b+','Markersize',10,'parent',handles.RMSaxes);
    plot(rmst(Jpnt),rrms(Jpnt),'m+','Markersize',10,'parent',handles.RMSaxes);

else
    plot(rmst,bsm(OPTIONS.selectLead,:),'k','parent',handles.RMSaxes);
    plot(rmst(beats),bsm(OPTIONS.selectLead,beats),'or','Markersize',10,'parent',handles.RMSaxes);
    
    plot(rmst(pwaves),bsm(OPTIONS.selectLead,pwaves),'kx','Markersize',10,'parent',handles.RMSaxes);
    plot(rmst(twaves),bsm(OPTIONS.selectLead,twaves),'b+','Markersize',10,'parent',handles.RMSaxes);
    plot(rmst(Jpnt),bsm(OPTIONS.selectLead,Jpnt),'m+','Markersize',10,'parent',handles.RMSaxes);
    
    minA = -max(abs(bsm(OPTIONS.selectLead,:)));
    maxA = -minA;
end
axis([rmst(1) rmst(end) minA maxA]);

% Overlay RMS with signal
os = get(handles.OverlaySignal,'Value');
siglist = get(handles.OverlaySignal,'String');
if os>1 % 1 is none
    sig=str2num(siglist{os});
    osig=DATA.ORG(sig,:);
    osig=osig- ( lowpassma(osig, 1 / ( 0.1 * DATA.sampT ) ) );
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
steps = get(handles.sliderSignal,'sliderstep');
OPTIONS.viewt0 = round( max(1, get(handles.sliderSignal,'Value') * size(DATA.BSM,2) ) );
OPTIONS.viewt0 = min(OPTIONS.viewt0, size(DATA.ORG,2) - round(OPTIONS.viewt0 + (steps(2) * size(DATA.ORG,2))));

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

% set(handles.plotAxes,'Position',plotp);
% set(handles.RMSaxes,'Position',rmsp);
% set(handles.uipanelControls,'Position',contp);


%% --- Executes on button press in pushZoomIn.
function pushZoomIn_Callback(hObject, eventdata, handles)
% hObject    handle to pushZoomIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global OPTIONS
global DATA
OPTIONS.viewt0 = round( max(1, get(handles.sliderSignal,'Value') * size(DATA.BSM,2) ) );
steps = get(handles.sliderSignal,'sliderstep');
steps = steps / 2;
set(handles.sliderSignal,'sliderstep',steps);

% OPTIONS.viewt0 = round(get(handles.sliderSignal,'Value') + (steps(2)/2 * size(DATA.ORG,2))); % PO
OPTIONS.viewt1 = min(size(DATA.ORG,2),round(OPTIONS.viewt0 + (steps(2) * size(DATA.ORG,2))));

plotBSM(handles)

%% --- Executes on button press in pushZoomOut.
function pushZoomOut_Callback(hObject, eventdata, handles)
% hObject    handle to pushZoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global DATA
global OPTIONS

OPTIONS.viewt0 = round( max(1, get(handles.sliderSignal,'Value') * size(DATA.BSM,2) ) );
steps = get(handles.sliderSignal,'sliderstep');
steps = steps * 2;
set(handles.sliderSignal,'sliderstep',steps);
% OPTIONS.viewt0 =max(0,round(get(handles.sliderSignal,'Value')-steps(2)/2*size(DATA.ORG,2)));%PO
OPTIONS.viewt1 = min(size(DATA.ORG,2),round(OPTIONS.viewt0 + (steps(2) * size(DATA.ORG,2))));
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
selMarker = find(abs(DATA.BEATS-selTime) < 50);
maxA=max(DATA.RMS);
if strcmp(get(handles.selectBSMFig,'SelectionType'),'extend')% add zeropoint
    if isempty(find(abs(DATA.BEATS-selTime) < 100)) % no marker found within 100 ms
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
    doBaseline(handles,0);
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
    plotBSM(handles)
elseif ~isempty(selMarker) && (OPTIONS.viewt1 - OPTIONS.viewt0) < 30000
    figpoint1=get(handles.selectBSMFig,'currentpoint');
    rect=[figpoint1(1,1), figpoint1(1,2), 2, maxA];
    rbbox;
    point1=get(hObject,'currentpoint');
    DATA.BEATS(selMarker) = round(point1(1)/DATA.sampT);
    doBaseline(handles,0);
else
    OPTIONS.selectLead = -1;
    if (OPTIONS.viewt1 - OPTIONS.viewt0) > 30000
        rbbox;
        point2 = get(hObject,'CurrentPoint');% button up detected
        xmax=max(point1(1,1),point2(1,1));
        xmin=min(point1(1,1),point2(1,1));
        if xmax < size(DATA.BSM,2) - 1 %add one second min width
            xmax = max(xmin+1,xmax);
        else
            xmin = min(xmin,xmax-1);
        end
        OPTIONS.viewt0 = max(1,round(xmin/DATA.sampT));
        OPTIONS.viewt1 = min(size(DATA.BSM,2),round(xmax/DATA.sampT));
        set(handles.sliderSignal,'Value',(OPTIONS.viewt0)/size(DATA.BSM,2));
        step = (OPTIONS.viewt1 - OPTIONS.viewt0)/size(DATA.BSM,2);
        steps = [step/10 step];
        set(handles.sliderSignal,'sliderstep',steps);
    end
    plotBSM(handles)   
end




% --- Executes on mouse press over axes background.
function plotAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to plotAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DATA
global OPTIONS
% 	yshift(j)=2*raster(2)+1-2*DATA.LAY(i,2);
% 	xshift(j)=DATA.LAY(i,1)-1;

point1=get(hObject,'currentpoint');
x=floor(1+point1(1,1));
y=-round((point1(1,2) - DATA.LAY(1,2) * 2-1)/2);
% disp(num2str([point1(1,1) point1(1,2) x y]))
signalIndex=find(DATA.LAY(2:end,1)==x & DATA.LAY(2:end,2)==y);

if strcmp(get(handles.selectBSMFig,'SelectionType'),'alt')%
    if ~isempty(signalIndex)
        OPTIONS.selectLead = signalIndex;
    else
        OPTIONS.selectLead = -1;
    end
    plotBSM(handles);
elseif ~isempty(signalIndex)
    OPTIONS.selectLead = -1;
    DATA.remove(signalIndex)=xor(DATA.remove(signalIndex),1);
    plotBSM(handles);
%     doBaseline(handles,0);   
end

% --- Executes on button press in checkboxHighPass.
function checkboxHighPass_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxHighPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxHighPass

% global DATA

% DATA.BSM=DATA.ORG;
% if get(hObject,'Value')
% 	DATA.BSM=DATA.BSM-lowpassma(DATA.BSM,round(DATA.len20ms/DATA.sampT)); % 5Hz
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
doBaseline(handles,1);








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
    bsm = DATA.BSMfilt;
else
    bsm = DATA.BSMraw;
end
ShowBSMonThorax(bsm(:,OPTIONS.viewt0:OPTIONS.viewt1),OPTIONS.leadsys); % plot zeromean data on thorax


% --- Executes on button press in readmatbutton.
function readmatbutton_Callback(hObject, eventdata, handles)
% hObject    handle to readmatbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

return
global DATA
global OPTIONS
global D mlas
[OPTIONS.FileName,OPTIONS.PathName,FilterIndex] = uigetfile({'*.mat','*.*'},'Select data file.','D:/Users/peter/');
disp(OPTIONS.FileName);

load(fullfile(OPTIONS.PathName,OPTIONS.FileName),'DATA','D','mlas');

set(handles.PQint,'String',num2str(DATA.PQint*DATA.sampT*1000));
OPTIONS.viewt0=1;
OPTIONS.viewt1=ceil(max(100,size(DATA.ORG,2)/4));
set(handles.sliderSignal,'sliderStep',[0.01,0.1])
set(handles.sliderSignal,'min',0);


OPTIONS.zeromean=0;
OPTIONS.zmchannels=[];
set(handles.checkboxHighPass,'value',0);
set(handles.sliderSignal,'value',0);

set(handles.checkboxZeromean,'value',OPTIONS.zeromean);
set(handles.zmchannels,'String','');

%Init OverlaySignal list
siglist=[{'(none)'};cellstr(num2str((1:size(DATA.ORG,1))'))];
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

