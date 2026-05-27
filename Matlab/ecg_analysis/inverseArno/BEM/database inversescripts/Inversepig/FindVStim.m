sub=9;

switch sub
    case 7
        egmrange=65:68;
        iaegm=0;
        iendoegm=3;
        iepiegm=1;
    case 8
        egmrange=65:68;
        iaegm=0;
        iendoegm=2; % could be 1
        iepiegm=4;
        
    case 9
        egmrange=65:76;
        iaegm=12;
        iendoegm=11;
        iepiegm=1;
    case 10
        egmrange=65:79;
        iaegm=9;
        iendoegm=1;
        iepiegm=5;
end


if ~exist('FileName','var') || isempty(FileName)
    [FileName,PathName,FilterIndex] = uigetfile('*.ecgspecs','Select files'...
        ,'/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats','MultiSelect','on');
end
if ~iscell(FileName)
    FileName={FileName};
end
for ifile=1:length(FileName);
    FilePath=(fullfile(PathName,FileName{ifile}));
    
    
    
    
    

end
