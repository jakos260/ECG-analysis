function ECG = readECG(filename,type,t0,t1)


% volume 1 Comprehensive electrocardiography pag. 384
% III = II - I;
% aVr= -0.5 * (I +II)       Vr = -1/3(I + II)
% aVl = I - 0.5 II          Vl = 2/3 * I - 1/3 II
% aVf = II - 0.5 I          Vf = 2/3 II - 1/3 I


% III = II - I;
% aVr = -1/2 * ( I + II );
% aVl = I - 0.5 * II;
% aVf = II - 0.5 * I;


if strcmpi(type,'mortara')
    ECG = readMortaraHere(filename,t0,t1);
elseif strcmpi(type,'heartsciences')
    ECG = readHeartSciences(filename);
elseif strcmpi(lower(type),'powermed')
    ECG = readPowermed(filename);
elseif strcmpi(type,'bin')
    ECG = readPrucka(filename,t0,t1);
elseif strcmpi(lower(type),'hfecg')
    ECG = readHfECG(filename);
elseif strcmpi(type,'erasmus')
    ECG = readErasmus(filename,t0,t1);   
elseif strcmpi(type,'larisa')
    ECG = readLarisa(filename);
elseif strcmpi(type,'lumc')
    ECG = readLUMC(filename);
elseif strcmpi(type,'txt')
    ECG = readTXT(filename);
elseif strcmpi(type,'ensite')
    ECG = readEnsite(filename);
elseif strcmpi(type,'cardiopcsv')
    ECG = readCardioPerfectCSV(filename);         
elseif strcmpi(type,'mortaracsv')
%     ECG = readMortaraCSV(filename);
elseif strcmpi(type,'viccsv')
    ECG = readVicCSV(filename);
elseif strcmpi(type,'fysiologic')
    ECG = readFysiologic(filename);
elseif strcmpi(type,'csvfysiologic')
    ECG = readLeadsCSV( filename);
elseif strcmpi(type,'elicsv')
    ECG = table2array(readtable(filename))';
    ECG = ECG(2:end,:) / 937.5;    
elseif strcmpi(type,'haga')
    if ~isempty(strfind(filename,'.bdf'))
        filename=filename(1:end-4);
    end
    H=ImportBDFHeader(filename);
    realname = [filename, '.bdf'];
    sizeHeader = 256 + 256 * H.channels;
    fin = fopen(realname,'r');
    fseek(fin,0,'eof');
    endOfFile = ftell(fin);
    fclose(fin);
    Exg1 = 66;%H.channels - 8;
    Exg2 = 65;%H.channels - 7;
    Exg3 = 67;%H.channels - 6;
    channels= [ 1:48 60 61 62 49 : 59 Exg1 Exg2 Exg3];
    chnim65=[10:64 3:9 1 2 65 ];
    if H.channels < 65
        error('not enough channels')
    end
    channel=1;
    sig = H.sensor.gain(channels(channel)) * ChannelReaderBDF(realname,H.channels,H.nSamples,H.nTrials,channels(channel),H.sampleRate,endOfFile,sizeHeader)';
    sig = sig(max(1,t0*H.sampleRate):t1*H.sampleRate);
    sig = resampleAvoHere(sig,1,length(sig),length(sig)*1000/H.sampleRate);
    SIGS=zeros(length(chnim65),length(sig));
    ichan = 1;
    SIGS(chnim65(ichan),:)=sig;
    % sensors zijn  7 kanalen
    % 64 kanalen waarvan kanaal 63 en 64 niet gebruikt
    % van de EXECG de eerste 3 kanalen
    hw=waitbar(0,'reading data');
    for channel=2:size(SIGS,1)%H.channels
        sig = H.sensor.gain(channels(channel)) * ChannelReaderBDF(realname,H.channels,H.nSamples,H.nTrials,channels(channel),H.sampleRate,endOfFile,sizeHeader)';
        sig = sig(max(1,t0*H.sampleRate):t1*H.sampleRate);
        SIGS(chnim65(channel),:) = resampleAvoHere(sig,1,length(sig),length(sig)*1000/H.sampleRate);
        waitbar(channel/size(SIGS,1),hw);drawnow
    end
    close(hw)
    ECG = SIGS / 1000;
elseif strcmpi(type,'bdf')
    D = readBDF(filename,[],t0,t1);
    nsamp = 1000 * (size(D.signals,1)) / D.Fs;
    x=1:size(D.signals,1);
    str=(nsamp-1)/(size(D.signals,1));
    x= 1 + (x-1)*str;
    ECG=interp1(x,D.signals,1 : nsamp,'cubic')';
elseif strcmpi(type,'mat')
    ECG = loadmat(filename);
else
    error('type not suported')
end
%%
function ECG = readdossantos(dirname)

fid=fopen(fullfile(dirname,'D1'),'r','b');I=fread(fid,inf,'float32');fclose(fid);
fid=fopen(fullfile(dirname,'D2'),'r','b');II=fread(fid,inf,'float32');fclose(fid);
fid=fopen(fullfile(dirname,'V1'),'r','b');V1=fread(fid,inf,'float32');fclose(fid);
fid=fopen(fullfile(dirname,'V2'),'r','b');V2=fread(fid,inf,'float32');fclose(fid);
fid=fopen(fullfile(dirname,'V3'),'r','b');V3=fread(fid,inf,'float32');fclose(fid);
fid=fopen(fullfile(dirname,'V4'),'r','b');V4=fread(fid,inf,'float32');fclose(fid);
fid=fopen(fullfile(dirname,'V5'),'r','b');V5=fread(fid,inf,'float32');fclose(fid);
fid=fopen(fullfile(dirname,'V6'),'r','b');V6=fread(fid,inf,'float32');fclose(fid);

aVr = -1/2 * ( I + II );
aVl = I - 0.5 * II;
aVf = II - 0.5 * I;
ECG = [I II (II-I) aVr aVl aVf V1 V2 V3 V4 V5 V6 ]';

%%

function ECG = readHeartSciences(filename)
A=xml2struct(filename);
leads=A.HeartTestLab.ECG_Export_Request.ECG_Result_Waveforms.lead;
ECG=[];
for i=1:length(leads)
    ecg = str2num(strrep(leads{i}.ecgLeadValueDigits.Text,',',' '));
    
    if strcmp(leads{i}.ecgLeadType.Text,'I')
        ECG(1,:) = ecg;
    elseif strcmp(leads{i}.ecgLeadType.Text,'II')
        ECG(2,:) = ecg;
    elseif strcmp(leads{i}.ecgLeadType.Text,'III')
        ECG(3,:) = ecg;
    elseif strcmp(leads{i}.ecgLeadType.Text,'aVR')
        ECG(4,:) = ecg;
    elseif strcmp(leads{i}.ecgLeadType.Text,'aVL')
        ECG(5,:) = ecg;
    elseif strcmp(leads{i}.ecgLeadType.Text,'aVF')
        ECG(6,:) = ecg;
    elseif strcmp(leads{i}.ecgLeadType.Text,'V1')
        ECG(7,:) = ecg;
    elseif strcmp(leads{i}.ecgLeadType.Text,'V2')
        ECG(8,:) = ecg;
    elseif strcmp(leads{i}.ecgLeadType.Text,'V3')
        ECG(9,:) = ecg;
    elseif strcmp(leads{i}.ecgLeadType.Text,'V4')
        ECG(10,:) = ecg;
    elseif strcmp(leads{i}.ecgLeadType.Text,'V5')
        ECG(11,:) = ecg;
    elseif strcmp(leads{i}.ecgLeadType.Text,'V6')
        ECG(12,:) = ecg;
    else
        error('stop');
    end
end

%%*************************************************************************
function ECG = readPowermed(filename)

txt = fileread(filename);
D=jsondecode(txt);
if isfield( D,'leads') && ...
    isfield(D.leads,'I') && isfield(D.leads,'II') && isfield(D.leads,'III') && ...
    isfield(D.leads,'aVR') && isfield(D.leads,'aVL') && isfield(D.leads,'aVF') && ...
    isfield(D.leads,'V1') && isfield(D.leads,'V2') && isfield(D.leads,'V3') && ...
    isfield(D.leads,'V4') && isfield(D.leads,'V5') && isfield(D.leads,'V6')    
    tmin= min([ length(D.leads.I.ecg),...
                length(D.leads.II.ecg),...
                length(D.leads.III.ecg),...
                length(D.leads.aVR.ecg),...
                length(D.leads.aVL.ecg),...
                length(D.leads.aVF.ecg),...
                length(D.leads.V1.ecg),...
                length(D.leads.V2.ecg),...
                length(D.leads.V3.ecg),...
                length(D.leads.V4.ecg),...
                length(D.leads.V5.ecg),...
                length(D.leads.V6.ecg)]);
    I = D.leads.I.ecg(1:tmin)';
    II= D.leads.II.ecg(1:tmin)';
    III= D.leads.III.ecg(1:tmin)';
    aVr = -1/2 * ( I + II );
    aVl = I - 0.5 * II;
    aVf = II - 0.5 * I;
    ecg = [ I;II;III;aVr;aVl;aVf;...
            D.leads.V1.ecg(1:tmin)';...
            D.leads.V2.ecg(1:tmin)';...
            D.leads.V3.ecg(1:tmin)';...
            D.leads.V4.ecg(1:tmin)';...
            D.leads.V5.ecg(1:tmin)';...
            D.leads.V6.ecg(1:tmin)'];
    R= rms(ecg); 
    if R(1) > 2 || any(isnan(R))
        ecg(:,isnan(R))=[];
    end
    fs = D.leads.I.metadata.fs;
    ECG = resampleAvo(ecg,1,size(ecg,2),round(size(ecg,2)*1000/320));
else
    error('not right format')
end
%%*************************************************************************
function ECG = readFysiologic(filename)

fid=fopen(filename,'r');
DATA = fread(fid,Inf,'uint8');
fclose(fid);

PACEPULSE=255;
ECG = zeros(9,10000);
t=ones(9,1);

Ch=1;
for i=535:length(DATA)
    Val=DATA(i); 
    if (Val < PACEPULSE ) %/ Val < 252
        if     (Ch==1) ECG(2,t(2))=255-Val;t(2)=t(2)+1; Ch=Ch+1;
        elseif (Ch==2) ECG(1,t(1))=Val;t(1)=t(1)+1;     Ch=Ch+1;
        elseif (Ch==3) ECG(3,t(3))=Val;t(3)=t(3)+1;     Ch=Ch+1;
        elseif (Ch==4) ECG(4,t(4))=255-Val;t(4)=t(4)+1; Ch=Ch+1;
        elseif (Ch==5) ECG(5,t(5))=255-Val;t(5)=t(5)+1; Ch=Ch+1;
        elseif (Ch==6) ECG(6,t(6))=255-Val;t(6)=t(6)+1; Ch=Ch+1;
        elseif (Ch==7) ECG(7,t(7))=255-Val;t(7)=t(7)+1; Ch=Ch+1;
        elseif (Ch==8) ECG(8,t(8))=255-Val;t(8)=t(8)+1; Ch=Ch+1;
        elseif (Ch==9) ECG(9,t(9))=255-Val;t(9)=t(9)+1; Ch=Ch+1;
        end
        if Ch>9
            Ch=1;
        end
    end

end
ECG=ECG(:,1:min(t)-1);

ECG(mean(ECG,2)==255,:)=0;

ECGh = ECG - lowpassma(ECG,100);
snr(ECGh,ECGh-lowpassma(ECGh,5))
ecg = ECGh - lowpassma(ECGh,5);
SNR = 10*log10(sum(ECGh(:).^2)/sum(ecg(:).^2));

aVr = -1/2 * ( ECG(1,:) + ECG(2,:) );
aVl = ECG(1,:) - 0.5 * ECG(2,:);
aVf = ECG(2,:) - 0.5 * ECG(1,:);

ECG = [ECG(1:3,:);[aVr; aVl; aVf] ;ECG(4:end,:); ];
ECG = ECG * 3.5/255; % We need to know the value per bit


% ecg = lowpassma(ECG,4);

ECG = resampleAvo(ECG,1,size(ECG,2),size(ECG,2)*5);
%%*************************************************************************
function ECG = readHfECG(filename,varargin)

%SignalPlant log file for:..
% Sampling frequency [Hz]:5000
% Datatype: 4-byte float
% Channels:57
% Samples:3000000
% Channel names and units:
a=strfind(filename,'.log');
fnBin = [filename(1:a(end)) 'bin'];
fidInf=fopen(filename,'rt');
fgetl(fidInf);
txt = fgetl(fidInf);
a= strfind(txt,'Sampling frequency [Hz]:');
sampleFreq = str2num(txt(a(end)+24:end));
fgetl(fidInf);
txt = fgetl(fidInf);
a= strfind(txt,'Channels:');
nrChannels = str2num(txt(a(end)+9:end));
txt = fgetl(fidInf);
a= strfind(txt,'Samples:');
nrSamples = str2num(txt(a(end)+8:end));
fgetl(fidInf);
index = zeros(12,1);
for i=1:14
    txt = fgetl(fidInf);
    A=strsplit(txt);    
    if strcmp(A{1},'I')
        index(1)=i;
    elseif strcmp(A{1},'II')
        index(2)=i;
    elseif strcmp(A{1},'III')
        index(3)=i;
    elseif strcmp(A{1},'aVR')
        index(4)=i;
    elseif strcmp(A{1},'aVL')
        index(5)=i;
    elseif strcmp(A{1},'aVF')
        index(6)=i;
    elseif strcmp(A{1},'V1')
        index(7)=i;
    elseif strcmp(A{1},'V2')
        index(8)=i;
    elseif strcmp(A{1},'V3')
        index(9)=i;
    elseif strcmp(A{1},'V4')
        index(10)=i;
    elseif strcmp(A{1},'V5')
        index(11)=i;
    elseif strcmp(A{1},'V6')
        index(12)=i;
    end    
end

fclose(fidInf);

tbegin = 1;
if nargin >= 2
    tbegin = max(max(1,t0), round(t0 * sampleFreq) );
end
tend = inf;
if nargin == 3
    tend = min(nrSamples,round(t1*sampleFreq));
end

fidbin=fopen(fnBin,'r');
if fidbin>=0   
    ECG= fread(fidbin,inf,'float32');
    fclose(fidbin);
    if length(ECG) == nrChannels * nrSamples
        ECG = reshape(ECG,nrChannels,length(ECG)/nrChannels);
        if sampleFreq > 1000
            ECG = lowpassma(ECG(index,:), round(sampleFreq / 1000)); % in uV
        end
        tbegin = min(size(ECG,2),tbegin);
        tend = min(size(ECG,2),tend);
        ECG = ECG(:,tbegin:tend)/1000;
        ECG = resampleAvoHere(ECG,tbegin,tend,(tend-tbegin+1)*1000/sampleFreq);
    else
        error(['could not read data ' fn]);
    end
end


%%*************************************************************************
function ECG = readPrucka(varargin)

fn = varargin{1};

% Patient =
% Description =
% Date = 1/10/2013 10:41:57
% Number of Channel = 12
% Points for Each Channel = 1468749
% Data Sampling Rate = 977 points/second
% Start Time = 1/10/2013 9:11:38
% Stop Time = 1/10/2013 9:36:42
% Units: mmHg for pressure and mV for all others
% Channel Number  Channel Label

a=strfind(fn,'.inf');
fnBin = [fn(1:a(end)) 'bin'];
fnTxt = [fn(1:a(end)) 'txt'];


fidInf=fopen(fn,'rt');
fgetl(fidInf);
fgetl(fidInf);
fgetl(fidInf);
txt = fgetl(fidInf);
a= strfind(txt,'=');
nrChannels = str2num(txt(a(end)+1:end));
txt = fgetl(fidInf);
a= strfind(txt,'=');
nrSamples = str2num(txt(a(end)+1:end));
txt = fgetl(fidInf);
a= strfind(txt,'=');
b= strfind(txt,'points/second');
sampleFreq = str2num(txt(a(end)+1:b-1));
fclose(fidInf);
tbegin = 1;
if nargin >= 2
    tbegin = max(1,round(varargin{2}*sampleFreq));
end
tend = inf;
if nargin == 3
    tend = round(varargin{3}*sampleFreq);
end

fidbin=fopen(fnBin,'r');
if fidbin>=0
    
    ECG= fread(fidbin,inf,'double');

    fclose(fidbin);
    if length(ECG) == nrChannels * nrSamples
        ECG=reshape(ECG,nrChannels,length(ECG)/nrChannels);
    %     ECG=ECG(4:end,:);
    %     ECG(1:3,:) = ECG(1:3,:)/1.5;
    %     
        tbegin = min(size(ECG,2),tbegin);
        tend = min(size(ECG,2),tend);
        ECG = resampleAvoHere(ECG,tbegin,tend,1000*(tend-tbegin+1)/sampleFreq);
    else
        error(['could not read data ' fn]);
    end
else 
    T = readtable(fnTxt);
    if isempty(T)
        error(['could not read data ' fn]);
    else
        ECG = T{:,1:12}';
        tbegin = min(size(ECG,2),tbegin);
        tend = min(size(ECG,2),tend);
        ECG = resampleAvoHere(ECG,tbegin,tend,1000*(tend-tbegin+1)/sampleFreq);
    end
end



%%**********************************************************************
function ECG = readErasmus(filename,t0,t1)
if strfind(filename,'.bwv')
    fid=fopen(filename);
else
    fid=fopen([ filename '.bwv']);
end
A=fread(fid,Inf,'int16');
fclose(fid);

fidInf=fopen([ filename '.ini']);
fgetl(fidInf); % patient id
fgetl(fidInf); % modality
fgetl(fidInf); % first sample
fgetl(fidInf); % lastsampl
txt = fgetl(fidInf);
a= strfind(txt,'=');
sampleFreq = str2num(txt(a(end)+1:end));
fgetl(fidInf); % duration
txt = fgetl(fidInf);
a= strfind(txt,'=');

nrChannels = str2num(txt(a(end)+1:end));
A = reshape(A,nrChannels,length(A)/nrChannels);
ECG= zeros(9,size(A,2));
keep = zeros(nrChannels,1);
for i=1:nrChannels
    fgetl(fidInf); % empty line
    fgetl(fidInf); % channleline
    label = fgetl(fidInf); % label
    fgetl(fidInf); % bytes per sample
    txt = fgetl(fidInf);
    a= strfind(txt,'=');
    sensitivity = str2num(txt(a(end)+1:end));
    A(i,:) = A(i,:) * sensitivity;
    fgetl(fidInf); % units
    fgetl(fidInf); % correction factor
    fgetl(fidInf); % baseline
    fgetl(fidInf); % timeskew
    
    if ~isempty(strfind(label,'aVR'))
        ECG(1,:)= A(i,:);
    elseif ~isempty(strfind(label,'aVL'))
        ECG(2,:)= A(i,:);
    elseif ~isempty(strfind(label,'aVF'))
        ECG(3,:)= A(i,:);
    elseif ~isempty(strfind(label,'V1'))
        ECG(4,:)= A(i,:);
    elseif ~isempty(strfind(label,'V2'))
        ECG(5,:)= A(i,:);
    elseif ~isempty(strfind(label,'V3'))
        ECG(6,:)= A(i,:);
    elseif ~isempty(strfind(label,'V4'))
        ECG(7,:)= A(i,:);
    elseif ~isempty(strfind(label,'V5'))
        ECG(8,:)= A(i,:);
    elseif ~isempty(strfind(label,'V6'))
        ECG(9,:)= A(i,:);
    end
end
ECG(1:3,:)= ECG(1:3,:)/1.5;
ECG = zeromean(ECG);
ECG = resampleAvoHere(ECG,1,size(ECG,2),1000*(size(ECG,2))/sampleFreq);

%%*************************************************************************
function ECG = readMortaraHere(varargin)

fn = varargin{1};
tbegin = 1;
if nargin >= 2
    tbegin = max(1,round(varargin{2}*1000.0));
end
tend = inf;
if nargin == 3
    tend = round(varargin{3}*1000);
end
% Add the missing Lead information
% aVl = I - 0.5 II
% aVr= -0.5 * (I +II)
% aVf = II - 0.5 I

fid=fopen(fn,'r');
ECG=fread(fid,'int16')*15/16000;
fclose(fid);
ECG=reshape(ECG,8,length(ECG)/8);
tend = min(tend,size(ECG,2));
ECG=ECG(:,tbegin:tend);

aVr = -1/2 * ( ECG(1,:) + ECG(2,:) );
aVl = ECG(1,:) - 0.5 * ECG(2,:);
aVf = ECG(2,:) - 0.5 * ECG(1,:);

ECG = [ECG(1,:);ECG(2,:); ECG(2,:)-ECG(1,:);[aVr; aVl; aVf] ;ECG(3:end,:); ];
%%
function ECG = readLarisa(fn)
fid=fopen(fn,'rt');
a=fgetl(fid);
a=fgetl(fid);
a=fgetl(fid);
a=fgetl(fid);
a=fgetl(fid);
a=fgetl(fid);
ECG=zeros(5000,12);
i=0;
while ~feof(fid)
    i=i+1;
    str=fgetl(fid);
    str=strrep(str,',',' ');
    ECG(i,:)=str2num(str);
end
fclose(fid);
ECG=ECG'/1000;
ECG = resampleAvoHere(ECG,1,length(ECG),length(ECG)*2);

%% -------------------- readLUMC ---------------------------------------
function ECG = readLUMC(fn)
fid=fopen(fn,'rt');
a=fgetl(fid);
ECG=zeros(100000,8);
i=0;
while ~feof(fid)
    i=i+1;
    str=fgetl(fid);
    str=strrep(str,',',' ');
    ECG(i,:)=str2num(str);
end
fclose(fid);
ECG=ECG'/1000;
ECG= ECG(:,1:i);

aVr = -1/2 * ( ECG(1,:) + ECG(2,:) );
aVl = ECG(1,:) - 0.5 * ECG(2,:);
aVf = ECG(2,:) - 0.5 * ECG(1,:);
ECG = [ECG(1,:);ECG(2,:); ECG(2,:)-ECG(1,:);[aVr; aVl; aVf] ;ECG(3:end,:); ];
ECG = resampleAvoHere(ECG,1,length(ECG),length(ECG)*2);

%% ---------------------------- readTXT --------------------------------
function ECG = readTXT(filename)
fid=fopen(filename,'rt');
a = length(str2num(fgetl(fid)));
fclose(fid);
fid=fopen(filename,'rt');
ECG = cell2mat(textscan(fid,'%f'));
fclose(fid);
ECG=reshape(ECG,a,length(ECG)/a);
ECG(4:6) = ECG(4:6)/1.5;
ECG=ECG(4:12,:);

%% --------------------------readEnsite ----------------------------------
function ECG = readEnsite(fn)
fid=fopen(fn);
while  ~feof(fid)
    str =fgetl(fid);
    if contains(str,'Number of samples (rows)')
        str =fgetl(fid);
        [headers]=split(str,',');
        str = fscanf(fid,'%s');
        fclose(fid);
        useheaders = zeros(size(headers));
        for i=1:length(headers)
            if strcmp(headers{i},'I') ||...
               strcmp(headers{i},'II') ||...
               strcmp(headers{i},'III') ||...
               strcmp(headers{i},'aVR') ||...
               strcmp(headers{i},'aVL') ||...
               strcmp(headers{i},'aVF') ||...
               strcmp(headers{i},'V1') ||...
               strcmp(headers{i},'V2') ||...
               strcmp(headers{i},'V3') ||...
               strcmp(headers{i},'V4') ||...
               strcmp(headers{i},'V5') ||...
               strcmp(headers{i},'V6') 
                useheaders(i) = 1;
            elseif strcmp(headers{i},'t_ref') 
                refTimeCol = i;            
            end
        end        
        data =split(str,',');
        DATA = zeros(size(headers,1),(size(data,1)-1) / size(headers,1));                
        k=1;
        p=1;
        for i = 1 : size(DATA(:))
            if useheaders(p) || p == refTimeCol
                DATA(p,k)=str2double(data{i});
            end
            
            if rem(i,size(headers,1)) == 0
                k = k+1;
                p=1;
            else
                p=p+1;
            end            
        end
        ECG=DATA(useheaders==1,:);
        tmax = max(DATA(refTimeCol,:));
%         fs = mean((diff(DATA(refTimeCol,1:end))));
        ECG=resampleAvoHere(ECG,1,size(ECG,2),tmax*1000);
        break;
    end
end

%%
function PSI=resampleAvoHere(PHI,tbeg,tend,nsamp)
dim=size(PHI);
incr=nsamp/(tend-tbeg);
x=1:tend-tbeg+1;

str=(nsamp-1)/(tend-tbeg);
x=1+(x-1)*str;
xi=1:nsamp;
PSI=interp1(x,PHI(:,tbeg:tend)',xi,'pchip');
PSI=PSI';
%%
function     ECG12 = readVicCSV(fn)
% data is 1000Hz 1 uV
T = readtable(fn);
ECG = T{:,1:8}/1000;
use = zeros(size(ECG));
for i=2:size(T,1)
    str = T{i,9}{1};   
    use(i,:) = strcmp(strsplit(str(1:end-1),'.'),'T');    
end
ValidData = sum(use,2)>1;
ECG(use==0)=0;

ECG = ECG(ValidData,:)';
aVr = -1/2 * ( ECG(1,:) + ECG(2,:) );
aVl = ECG(1,:) - 0.5 * ECG(2,:);
aVf = ECG(1,:) - 0.5 * ECG(2,:);
ECG12 = [ECG(1,:);ECG(2,:); ECG(2,:)-ECG(1,:);[aVr; aVl; aVf] ;ECG(3:end,:); ];

%%
function     ECG12 = readCardioPerfectCSV(fn)

fid=fopen(fn);
str =fgetl(fid);
str =fgetl(fid);
headers=split(str,',');
sampleFreq = str2num(headers{4});
nvperLSB = str2num(headers{5})* 1e-6;

str =fgetl(fid);
str =fgetl(fid);
% headers=split(str,',');
fclose(fid);

ECG = csvread(fn,4,0)' * nvperLSB;

aVr = -1/2 * ( ECG(7,:) + ECG(8,:) );
aVl = ECG(7,:) - 0.5 * ECG(8,:);
aVf = ECG(8,:) - 0.5 * ECG(7,:);
ECG12 = [ECG(7,:);ECG(8,:); ECG(8,:)-ECG(7,:);[aVr; aVl; aVf] ;ECG(3:end,:); ];
%%

function     ECG12 = readMortaraCSV(fn)


ECG = csvread(fn)';

aVr = -1/2 * ( ECG(7,:) + ECG(8,:) );
aVl = ECG(7,:) - 0.5 * ECG(8,:);
aVf = ECG(8,:) - 0.5 * ECG(7,:);
ECG12 = [ECG(7,:);ECG(8,:); ECG(8,:)-ECG(7,:);[aVr; aVl; aVf] ;ECG(3:end,:); ];

function     ECG12 = readLeadsCSV(fn)
ECG=[];

fid=fopen(fn);
head = fgetl(fid);
while ~feof(fid) 
    str =fgetl(fid);
    ECG = [ ECG; str2num(sscanf(str,'%s'))];
end
fclose(fid);
ECG = ECG'/1000;

aVr = -1/2 * ( ECG(7,:) + ECG(8,:) );
aVl = ECG(7,:) - 0.5 * ECG(8,:);
aVf = ECG(8,:) - 0.5 * ECG(7,:);
ECG12 = [ECG(7,:);ECG(8,:); ECG(8,:)-ECG(7,:);[aVr; aVl; aVf] ;ECG(3:end,:); ];
ECG12=resampleAvoHere(ECG12,1,size(ECG12,2),size(ECG12,2)*2);



