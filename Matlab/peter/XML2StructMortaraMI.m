function DATA = XML2StructMortaraMI(filename) 

A=xml2struct(filename);


DATA.TYPICAL = A.ECG.TYPICAL_CYCLE.Attributes;

names=fieldnames(DATA.TYPICAL);

for i =1:length(fieldnames(DATA.TYPICAL))
    eval(['a = (DATA.TYPICAL.' [names{i}] ');'] )
    a =str2double(a); 
    
    eval(['DATA.TYPICAL.' [names{i}] '= a;'] )
end

if isempty(A.ECG.Attributes.VENT_RATE)
DATA.VENT_RATE=0;
else
DATA.VENT_RATE=str2num(A.ECG.Attributes.VENT_RATE);
end
if isempty(A.ECG.Attributes.AVERAGE_RR)
DATA.meanRR=0;
else
DATA.meanRR=str2num(A.ECG.Attributes.AVERAGE_RR);
end

% str= A.ECG.CHANNEL{1}.Attributes.DATA;



% B=sscanf(str,'%x','inf');


varargout{1} = DATA;
