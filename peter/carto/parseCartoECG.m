function [ECG, channels]= parseCartoECG(filename)

fid=fopen(filename);

str=fgetl(fid);
str=fgetl(fid);
if strfind(str,'Raw ECG to MV (gain) =')
    gain= str2num(str(length('Raw ECG to MV (gain) = '):end));
end
str=fgetl(fid);
str=fgetl(fid);
channels=[];
i=1;

parseString=[];
while ~isempty(str)
    [channel,str]=strtok(str,' ');
    if ~isempty(channel)
        channels{i} = channel;
        parseString=[parseString '%d ']; 
        i=i+1;
    end
end
CAR = textscan(fid,parseString );

ECG= double(cell2mat(CAR(1:end)))*gain;
ECG= ECG';

fclose(fid);