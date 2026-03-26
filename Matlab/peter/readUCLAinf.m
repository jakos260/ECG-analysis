function ECG = readUCLAinf(varargin)

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


