function ViewandExport(filepath,ElecPos,defsel,defstart,defend)
close all;
% clear all;
% clear global
%flags
fl_waterfall=0; % create waterfall plot
fl_sigplot=1; % plot signal
fl_selectbeats=1; % call Peters selectbeats
fl_onlyinmla=false;
if ~exist('ElecPos')
    ElecPos='pigbsm'; % % 'pbsm', 'amsbdf', 'amszip' , 'pigbsm', 'pigbsmegm': STW pig BSM met 4xEGM 'pigbsm09' '9lead' or ElecPos.
end
% filepath='/Users/peteroosterhoff/Documents/Werk/UCLA/proCIPS/proCIPS01/ECGcases/clinicalPVC1/ECG/X000B.cipsecg';
% filepath='/Users/peteroosterhoff/Documents/Werk/Brugada/Ajm009_5346982/BSPM/hagen/hagen2.bdf';
% filepath='/Users/peteroosterhoff/Documents/Werk/Brugada/Ajm009_5346982/BSPM/hagen/25-04-2005/hagen 25-04-2005 13_22_52.zip';
% filepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig03/BSM/varken2443.bdf';
%filepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig03/BasketEnSok/varken2443_3BasketEnSok+002.bdf'; %OSX
%filepath='Y:/Documents/Werk/STW/Data/Pig03/BasketEnSok/varken2443_3BasketEnSok+003.bdf'; % Windows
% filepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Bucket03/18grp1_EMMER_200312_VM+024.bdf';
% filepath='/Users/peteroosterhoff/Documents/Werk/AZU/LAMAP/Data/LAMAP04/ECG/LAMAp004.bdf';
% filepath='/Users/peteroosterhoff/Documents/Werk/AZU/LAMAP/Data/LAMAP05/ECG/LAMAP005NaCT.bdf';
% filepath='/Users/peteroosterhoff/Documents/Werk/AZU/LAMAP/Data/LAMAP02/ECG/LAMAP002.bdf';
if ~exist('filepath')
%     defaultpath='/Users/peteroosterhoff/Documents/Werk/STW/Klinische studie/AMC/Data/NICE001/Biosemi';
         defaultpath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig08/Biosemi/';
%     defaultpath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig08/Biosemi/';
%     defaultpath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig10/Biosemi/';
%     defaultpath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig10/Biosemi/';
% defaultpath='/Users/peteroosterhoff/Documents/Werk/UCLA/proCIPS/proCIPS01/ECGcases/clinicalPVC1/ECG';
%     defaultpath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Bucket01/';
%     defaultpasth='/Users/peteroosterhoff/Documents/Werk/STW/Data/Bucket03/';
%     defaultpath='/Users/peteroosterhoff/Documents/Werk/ARVC/DATA/ARVCp02/ECG/MUSE_20131211_155505_67000.csv';
%     defaultpath='/Users/peteroosterhoff/Documents/Werk/Brugada/DATA/Measurements/AJM060/BSPM/';
    [FileName,PathName,FilterIndex] = uigetfile({'*.bdf;*.zip;*.XML;*.xml;*.csv;*.cipsecg'},'Select data file.',defaultpath);
%     [FileName,PathName,FilterIndex] = uigetfile('*.bdf;*.zip','Select Biosemi data file.','/Volumes/Alzheimer/DATA_HD_Andre/Ajm048_6159323/BSPM/6159323/18-04-2006');
    filepath=fullfile(PathName,FileName);
end
% '/Users/peteroosterhoff/Documents/Werk/Brugada/DATA/Measurements/AJM059/BSPM');%

[folder, filename, ext]=fileparts(filepath);







% DefAns{1}='1:76';
defsel='1:78';
% chan=eval(DefAns{1});
if strcmpi(ext,'.bdf')
    if exist('defsel')
        DefAns{1}=defsel;
    else
        DefAns{1}='1:78';%'33:96'%'1:76';
    end
    if exist('defstart')
        DefAns{2}=defstart;
    else
        DefAns{2}='0'%'10';
    end
    if exist('defend')
        DefAns{3}=defend;
    else
        DefAns{3}='+60' %'+30';
    end
    if nargin>0
        D=Readbdf(filepath, DefAns{1}, DefAns{2},DefAns{3},'noui');
    else
        D=Readbdf(filepath, DefAns{1}, DefAns{2},DefAns{3});
    end
    D.filepath=filepath;
elseif strcmpi(ext,'.zip')
    if exist('defsel')
        DefAns{1}=defsel
    else
        DefAns{1}='1:65';
    end
    if exist('defstart')
        DefAns{2}=defstart;
    else
        DefAns{2}='13:23:02';
    end
    if exist('defend')
        DefAns{3}=defend;
    else
        DefAns{3}='+30';
    end
    dlg_title='Bdf-zip import dialog.';
    prompt{1}=sprintf('Select channels (format: 1:5 14). Default %s.',DefAns{1});
    prompt{2}='Start time hh:mm:ss ';
    
    prompt{3}='duration in seconds.';
    inputstr=inputdlg(prompt,dlg_title,1,DefAns);
    SelectChannels=eval(['[' inputstr{1} ']']);
    
    zipstarttime=inputstr{2};
    zipduration=inputstr{3};
    [D,mchd]=timeload(folder,zipstarttime,zipduration);
    %     D.Ts=mchd.RDAT.SampTime/1000; % us to ms ;
    D.signals=D.signals(:,SelectChannels);
    D.SelectChannels=SelectChannels;
elseif strcmpi(ext,'.xml')
    D=ReadMUSEXML(filepath);
    D.filepath=filepath;
    if length(D.chanlist)==8 && ~any(strcmpi('III',D.chanlist))
        warning('8 channels from XML.  Recreating III');
        D.chanlist={D.chanlist{1:2} 'III', D.chanlist{3:end}};
        D.signals(:,4:end+1)=D.signals(:,3:end);
        D.signals(:,3)=D.signals(:,2)-D.signals(:,1);
        warning('8 leads to 9WTC leads not implemented');
    end
     
elseif strcmpi(ext,'.csv');
    A=csvread(filepath,1,0)'/1000; % 12lead uV, 500Hz
    B=[A(4:6,:)/1.5 ; A(7:12,:)];
    D.signals = resample(B,1,size(B,2),size(B,2)*1/500*1000)';
    D.Ts=1;
    D.SelectChannels=1:9;
    D.chanlist={'VR','VL','VF','V1','V2','V3','V4','V6'};
    D.startpos=0;
    D.endpos=size(D.signals,2)/1000-1;
    D.firstsample=1;
elseif strcmpi(ext,'.cipsecg');
    A=loadmat(filepath); % 12lead 
    B=[A(4:6,:)/1.5 ; A(7:12,:)];
    D.signals = resample(B,1,size(B,2),size(B,2)*1/977*1000)';
    D.Ts=1;
    D.SelectChannels=1:9;
    D.chanlist={'VR','VL','VF','V1','V2','V3','V4','V6'};
    D.startpos=0;
    D.endpos=size(D.signals,2)/1000-1;
    D.firstsample=1;
else
    error('unknown file type');
end




%%
% Fs=2048;
% Ts=1000/Fs; % Sample time in ms
winstart=-200; % start of window in ms, negative is before event
winend=400;% end of window in ms, negative is before event
offsets=round(-200/D.Ts); % set signal to 0 x ms before event.
winstarts=floor(winstart/D.Ts);
winends=ceil(winend/D.Ts);


Fs=1000/D.Ts;
firstsample=Fs*D.startpos;
% lastsample=Fs*endtime;
lastsample=firstsample+size(D.signals,1)-1; % prevent rounding errors, accept roundig erro in startsample

if fl_waterfall
    load(fullfile(folder,'Events.mat')); %read markers into varaible 'events'
    firstevent=find((events>=firstsample-winstarts)&(events>=firstsample-offsets),1,'first');
    lastevent=find(events<=lastsample-winends,1,'last');
end

% [Data,MCHD]=ml_load(filepath,'SelectChannels',chan,'FirstSample',firstsample,'LastSample',lastsample); %
%[Data,MCHD]=ml_load(filepath,'SelectChannels',[26]);

% remove noise and humm






if fl_waterfall
    % make waterfall plot
    X=(winstarts:winends)*D.Ts;
    Y=events(firstevent:lastevent)/Fs/60; %time of events in minutes
    Z=zeros(length(Y),length(X));
    ii=1;
    for i=firstevent:lastevent
        Z(ii,1:length(X))=D.signals((winstarts:winends)+events(i)-firstsample,1)-D.signals(events(i)-firstsample+offsets,1);
        ii=ii+1;
    end
    % X=squeeze(X);
    % Y=squeeze(Y);
    figure;
    set(gcf,'Renderer','OpenGL');
    waterfall(X,Y,Z);
end

if fl_sigplot
    figure
    XX=D.startpos+(1:size(D.signals,1))*D.Ts/1000;
    xlim([D.startpos,D.startpos+10]);
    plot(XX,D.signals(:,2));
    drawnow
    figure % dummy figure
end

%%
if fl_selectbeats
    %global Nchannels;
    %Nchannels=60; % of toch 64?
    % psm; % maplab BSM layout
    
    
    [mla, standardc,leadsys]=elecpos2mlas( D.SelectChannels,ElecPos); % set ElecPos at start of this m-file.
    % mla seq L-R, T-D. org bdf chan numbers. standardc vr,vl,vf (Nijmegen), org bdf chan numbers.
    mlas=[]; % mla for selected channels only. First line with grid size added later
    delchan=[];
    standard=zeros(1,9);
    
    if fl_onlyinmla % remove signals not in mla layout
        for i=1:length(D.SelectChannels)
            pos=find(mla(:,3)==D.SelectChannels(i));
            if pos
                mlas(end+1,:)=[mla(pos,1:2),size(mlas,1)+1];
                % put in order of used selectchannels (=order of D.signals), numbers as after deleting unused signals
                poss=find(D.SelectChannels(i)==standardc,1,'first');
                if poss
                    standard(poss)=mlas(end,3);
                    % translate standard leads in org bdf channel numbers to
                    % signal in D.signals after deleting unused channels,
                    % vr,vl,vf order (Nijmegen)
                end
            else
                delchan(end+1)=i;
            end
        end
    else
        mlas=mla;
    end
    % Remove unreferenced channels
    if ~isempty(delchan)
        D.SelectChannels(delchan)=[];
        D.signals(:,delchan)=[];
    end
    
    if mlas(1,3)==0
        mlas=mlas(2:end,:); % remove grid size if present.
    end
    
    mlas(:,1)=mlas(:,1)-min(mlas(:,1))+1; % move grid to upper left
    mlas(:,2)=mlas(:,2)-min(mlas(:,2))+1;
    mlas=[max(mlas(:,1:2)),0; mlas]; % Now add first line with grid size
    if length(standard)==9 && all(standard>0)
        
        ind=1;% first line grid size
        for i=1:length(standard)
            pos=find(standard(i)==mlas(:,3),1,'first');
            if ~isempty(pos)  && ind<size(mlas,1)
                mlas=[mlas(1:ind,:);mlas(pos,:);mlas(ind+1:end,:)];
                mlas(pos+1,:)=[];
                ind=ind+1;
            end
        end
        
        % special case for amsterdam geo, 65 leads (and standard leads: else we would not be here)
        if max(mlas(:,3)==65) && strcmp(ElecPos(1:3),'ams')
            display('Assigning VL, VR and VF to channels fitting Ams geometry');
            dest=[64,63,65];% RV,LV, VF -> LV RV VF
            D.signals(:,dest)=D.signals(:,standard(7:9));
            mlas(8:10,3)=dest; %first line is grid size
            standard(7:9)=dest;
            
        end
    end
    
    % mlas=mla(ismember(mla(:,3),D.SelectChannels),:);
    % mlas(:,3)=1:size(mlas,1); % renumber channel numbers
    
    
    
    sig1000=interp1(1:size(D.signals,1),D.signals,1:1/D.Ts:size(D.signals,1),'PCHIP'); % D.Ts in ms, to Ts=0.001 s
    
    waitfor( selectBSMbeat(sig1000','lay',mlas,'sampt',1/1000,'leadsys',leadsys));
    global DATA
    response=inputdlg({'Short description (empty aborts save)','Comment'});
    if ~isempty(response) && ~isempty(response{1})
        D.comment=response{2};
        folderexp=[folder filesep 'export'];
        x=exist(folderexp,'file');
        if x==2
            folderexp=folder; % there is a file named 'export'
        elseif x~=7
            mkdir(folderexp);
        end
        save(fullfile(folderexp ...
            ,[filename '_' num2str(D.startpos) '_' num2str(D.endpos) '_'  response{1} '_' datestr(now,30) '.mat'])...
            ,'DATA','D','mlas');
    end
end



function [D, mchd]=timeload(folder,timstr,sz)
% adapted from bs_timeload

if nargin<3
    sz=10; % 10s default
end
tim=mod(datenum(timstr),1);
d=dir(fullfile(folder,'*.zip'));
for k=1:length(d)
    t(k)=mod(datenum(strrep(d(k).name(end-11:end-4),'_',':')),1);
end
fnr=find(t<=tim,1,'last');

if ~isempty(fnr)
    %get Sample interval. Read 1 sample
    [data,mchd]=ml_load(fullfile(folder,d(fnr).name),'FirstSample',1,'LastSample',1,...
        'MaxSamples',67108864,'Quite',0,'RangeType','Block','RecNr',1,'SelectChannels',[]);
    
    switch mchd.RDAT.SampTime
        case 488
            fs=2048;
        otherwise
            fs=round(1e6/mchd.RDAT.SampTime);
    end
    D.Ts=1000/fs;
    
    pos=1+floor((tim-t(fnr))*86400)*fs; % oostep1: 1+floor
    endpos=pos+round((str2double(sz)-1)*fs);
    [data,mchd]=ml_load(fullfile(folder,d(fnr).name),'FirstSample',pos,'LastSample',endpos,...
        'MaxSamples',67108864,'Quite',0,'RangeType','Block','RecNr',1,'SelectChannels',[]);
    D.signals=data';
    D.startpos=round(pos/fs);
    D.endpos=round(endpos/fs);
    D.filepath=fullfile(folder,d(fnr).name);
else
    error('No recording found at time %s in folder $s',timstr, folder);
end
