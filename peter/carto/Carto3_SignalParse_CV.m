function [DataPos, Electrograms, sr, EGMheads, numACQ, dataPath] = Carto3_SignalParse_CV(origPath)

% Author: Josh Blauer
% University of Utah
% Last Revision: Nov. 26, 2013
% Software parses the .txt files that are output from Carto3

data=[];
II=[];
CarPTS=[];

cd(origPath)

% Read in CAR.txt file for study of interest and get path to data directory.
[Car_filename,DataDir] = uigetfile('*car.txt','Select the car .txt file of interest (Carto3)');
dataPath = DataDir;
%Kick out of function if cancel is selected during uigetfile
if Car_filename==0;
    sr=0;
    DataPos = 0;
    Electrograms=0;
    sr=0;
    EGMheads=0;
    numACQ=0;
    return;
end
cd(origPath)

tic;

%Carto's sampling rate is 1000 Hz
sr=1000; %in Hz

% Read in *car.txt files
fid=fopen(strcat(DataDir,Car_filename));
CAR = textscan(fid, '%*c %*n %n %*n %n %n %n %n %n %n %n %n %n %*n %*n %n %*[^\n]','HeaderLines',1);
% Get point numbers
CarPTS=CAR{1,1};
% Get XYZ coords of CARTO points
CarCoords = [CAR{1,2} CAR{1,3} CAR{1,4}]; % These are the coords associate with the mapping device
fclose(fid);
[s,e] = regexp(Car_filename,'car');
% Look at first positions file to determine which catheter is being used to
% map, i.e. mapping and ablation catheter, lasso catheter, or penta-array
temp = dir([DataDir,Car_filename(1:s-1),'P',num2str(CarPTS(1)),'_','*Eleclectrode_Positions_OnAnnotation.txt']);
for k = 1:length(temp)
    [~,spntstr] = regexp(temp(k,1).name,[Car_filename(1:s-1),'P',num2str(CarPTS(1)),'_']);
    [~,epntstr] = regexp(temp(k,1).name,'CONNECTOR');
    pntstr = temp(k,1).name(spntstr+1:epntstr);
    DELIMITER = '\t';
    HEADERLINES = 2;
    mappingdeviceCoords = importdata([DataDir,Car_filename(1:s-1),'P',num2str(CarPTS(1)),....
        '_',pntstr,'_','Eleclectrode_Positions_OnAnnotation.txt'],DELIMITER, HEADERLINES);
    
    [~,DeviceNum] = intersect(mappingdeviceCoords.data(:,3:5),CarCoords(1,:),'rows');
    
    if ~isempty(dir([DataDir,Car_filename(1:s-1),'P',num2str(CarPTS(k)),....
        '_',pntstr,'_','Sensor_Positions_OnAnnotation.txt']))
        mappingsensorCoords = importdata([DataDir,Car_filename(1:s-1),'P',num2str(CarPTS(k)),....
            '_',pntstr,'_','Sensor_Positions_OnAnnotation.txt'],DELIMITER, HEADERLINES);
    
        [~,SensNum] = intersect(mappingsensorCoords.data(:,3:5),CarCoords(1,:),'rows');
    else
        SensNum = [];
    end
    
    if ~isempty(DeviceNum)
        break
    elseif ~isempty(SensNum)
        error('Positions based on sensor, not electrodes');
    elseif k == length(temp) && isempty(DeviceNum)
        error('Could not determine mapping device');
    end
    
    clear spntstr epntstr DELIMITER HEADERLINES
    clear mappingdeviceCoords
end

% Get Prime points (First of simultaneously acquired points)
% Use position files to determine which electrograms will be unique
clear temp
temp = dir([DataDir,Car_filename(1:s-1),'P*_',pntstr,'_Eleclectrode_Positions_OnAnnotation.txt']);

m = 1;
count = 1;
Ms = 1;
origFiles{count,1} = ...
    [DataDir,Car_filename(1:s-1),'P',num2str(CarPTS(m)),'_',pntstr,'_Eleclectrode_Positions.txt'];
for k = 1:numel(temp)
    [status,result] = system(['diff -q --brief ',...
        [DataDir,Car_filename(1:s-1),'P',num2str(CarPTS(m)),'_',pntstr,'_Eleclectrode_Positions_OnAnnotation.txt'],' ',...
        [DataDir,Car_filename(1:s-1),'P',num2str(CarPTS(k)),'_',pntstr,'_Eleclectrode_Positions_OnAnnotation.txt']]);
    if status == 1
        m = k;
        Ms = [Ms m];
        count = count+1;
        origFiles{count,1} = ...
            [DataDir,Car_filename(1:s-1),'P',num2str(CarPTS(k)),'_',pntstr,'_Eleclectrode_Positions.txt'];
    end
end
numACQ = length(Ms);
clear count m k status e result ans



if ~isempty(regexp(pntstr,'MAGNETIC_20_POLE_A_CONNECTOR', 'once'));
    EGMs = '20A_';   
elseif ~isempty(regexp(pntstr,'NAVISTAR_CONNECTOR', 'once'));
    EGMs = 'M';
elseif ~isempty(regexp(pntstr,'MAGNETIC_20_POLE_B_CONNECTOR', 'once'));
    EGMs = '20B_';
elseif ~isempty(regexp(pntstr,'CS_CONNECTOR', 'once'));
    EGMs = 'CS';
else
    error('Unrecognized Mapping Device')
end
clear temp

% Determine if unipolar and bipolar EGMs were acquired
clear DELIMITER HEADERLINES
% Read in Electrograms
DELIMITER = ' ';
%HEADERLINES = 3;
HEADERLINES = 4;
newData2 = importdata([DataDir,Car_filename(1:s-1),'P',num2str(CarPTS(Ms(1))),'_ECG_Export.txt'],...
    DELIMITER, HEADERLINES);

% Check for Unipolar EGMs
UniEGMs = [];
UniEGMs = strmatch([EGMs,'1('],newData2.colheaders);
if ~isempty(UniEGMs)
    EGMhere(1) = 1;
else
    EGMhere(1) = 0;
end

% Check for Bipolar EGMs
BipEGMs = [];
BipEGMs = strmatch([EGMs,'1-',EGMs,'2('],newData2.colheaders);
if ~isempty(BipEGMs)
    EGMhere(2) = 1;
else
    EGMhere(2) = 0;
end

% Check for limb lead II
IIcol = strmatch('II(',newData2.colheaders);

if sum(EGMhere) == 2
        EGMheads = newData2.colheaders([IIcol,UniEGMs:UniEGMs+9,BipEGMs:BipEGMs+8]); 
    elseif EGMhere(1) ==1 && EGMhere(2)==0
        EGMheads = newData2.colheaders([IIcol,UniEGMs:UniEGMs+9]);
    elseif EGMhere(1) ==0 && EGMhere(2)==1
        EGMheads = newData2.colheaders([IIcol,BipEGMs:BipEGMs+8]);
end

clear newData2 data colheaders

%% Find what device is mapping by looking for the electrode positions file.

%Parsing algorithm
wb=waitbar(0,'Please Wait...','Name','Loading Data');
for k=1:numACQ
    
    % Read in Electrode Positions (on annotation) file.
    
    DELIMITER = '\t';
    HEADERLINES = 2;
    newData = importdata(origFiles{k,1},DELIMITER, HEADERLINES);
    
    DataPos{k,1} = newData.data;
    clear DELIMITER HEADERLINES newData
    % Read in Electrograms
    DELIMITER = ' ';
    %HEADERLINES = 3;
    HEADERLINES = 4;
    newData2 = importdata([DataDir,Car_filename(1:s-1),'P',num2str(CarPTS(Ms(k))),'_ECG_Export.txt'],...
        DELIMITER, HEADERLINES);
    
    
    II=newData2.data(:,IIcol)*.003; %Gain *.003
    
    if EGMhere(1) ==1
    Uni=newData2.data(:,UniEGMs:UniEGMs+9)*.003; %Gain *.003
    end
    
    if EGMhere(2) ==1
    Bip=newData2.data(:,BipEGMs:BipEGMs+8)*.003; %Gain *.003
    end
    
    if sum(EGMhere) == 2
        Electrograms{k,1} = [II,Uni,Bip]; 
    elseif EGMhere(1) ==1 && EGMhere(2)==0
        Electrograms{k,1} = [II,Uni];
    elseif EGMhere(1) ==0 && EGMhere(2)==1
        Electrograms{k,1} = [II,Bip];
    end
    
    clear newData2
    waitbar(k/length(numACQ));
end
delete(wb);



t=toc;
end
