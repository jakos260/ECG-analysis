fprintf(['\n Starting ' mfilename '\n']);
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select export files for conversion to Nijmegen .mat'...
        ,'/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig10/BioSemi/export','MultiSelect','on');
if ~iscell(FileName)
    FileName={FileName};
end

for ifile=1:length(FileName)
   S=load(fullfile(PathName,FileName{ifile}));
   pathout=fullfile(PathName,'Nijmegenmat');
   if ~exist(pathout,'dir');
       mkdir(pathout);
   end
   savemat(fullfile(pathout,FileName{ifile}),S.DATA.ORG(1:64,:));
   
   
   
end;
fprintf(['\n Finished ' mfilename '\n']);