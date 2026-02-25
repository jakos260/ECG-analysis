function ECG = readECG(filename,type,t0,t1)


% volume 1 Comprehensive electrocardiography pag. 384
% aVr= -0.5 * (I +II) Vr = -1/3(I + II)
% aVl = I - 0.5 II    Vl = 2/3 * I - 1/3 II
% aVf = II - 0.5 I    Vf = 2/3 II - 1/3 I


if strcmpi(type,'mortara')
    ECG = readMortaraHere(filename,t0,t1);
elseif strcmpi(type,'bin')
    ECG = readPrucka(filename,t0,t1);
elseif strcmpi(type,'erasmus')
    ECG = readErasmus(filename,t0,t1);   
elseif strcmpi(type,'lumc')
    ECG = readLUMC(filename);
elseif strcmpi(type,'txt')
    ECG = readTXT(filename);
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
    sig = resample(sig,1,length(sig),length(sig)*1000/H.sampleRate);
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
        SIGS(chnim65(channel),:) = resample(sig,1,length(sig),length(sig)*1000/H.sampleRate);
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
fidbin=fopen(fnBin,'r');

tbegin = 1;
if nargin >= 2
    tbegin = max(1,round(varargin{2}*sampleFreq));
end
tend = inf;
if nargin == 3
    tend = round(varargin{3}*sampleFreq);
end

ECG= fread(fidbin,inf,'double');

fclose(fidbin);
if length(ECG) == nrChannels * nrSamples
    ECG=reshape(ECG,nrChannels,length(ECG)/nrChannels);
    ECG=ECG(4:end,:);
    ECG(1:3,:) = ECG(1:3,:)/1.5;
    
    tbegin = min(size(ECG,2),tbegin);
    tend = min(size(ECG,2),tend);
    ECG = resample(ECG,tbegin,tend,1000*(tend-tbegin+1)/sampleFreq);
else
    error(['could not read data ' fn]);
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
ECG = resample(ECG,1,size(ECG,2),1000*(size(ECG,2))/sampleFreq);

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

fid=fopen(fn,'r');
ECG=fread(fid,'int16')*15/16000;
fclose(fid);
ECG=reshape(ECG,8,length(ECG)/8);
tend = min(tend,size(ECG,2));
ECG=ECG(:,tbegin:tend);

Vr = -1/3 * ( ECG(1,:) + ECG(2,:) );
Vl = ECG(1,:) + Vr;
Vf = ECG(2,:) + Vr;
ECG1 = [ECG(1,:);ECG(2,:); Vf*1.5-Vl*1.5;[Vr; Vl; Vf] *1.5 ;ECG(3:end,:); ];
ECG = [Vr; Vl; Vf ;ECG(3:end,:); ];
ECG=ECG1;
%% -------------------- readLUMC ---------------------------------------
function ECG = readLUMC(fn)
fid=fopen(fn,'rt');
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
Vr = -1/3 * ( ECG(1,:) + ECG(2,:) );
% Vl =  2/3 * ECG(1,:) - 1/3 * ECG(2,:);
% Vf =  2/3 * ECG(2,:) - 1/3 * ECG(1,:);
Vl = ECG(1,:) + Vr; 
Vf = ECG(2,:) + Vr;
% ECG = [Vr;Vl;Vf; ECG(3:end,:); ECG(1,:); ECG(2,:) ];
ECG = [ECG(1:2,:); Vf-Vl; Vr; Vl; Vf; ECG(3:end,:)];
% ECG  = resampl(ECG,1,size(ECG,2),size(ECG,2)*2);

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
