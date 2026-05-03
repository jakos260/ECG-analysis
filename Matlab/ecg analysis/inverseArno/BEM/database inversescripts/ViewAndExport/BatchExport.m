clear
display(['Starting ' mfilename]);  
filename='/Users/peteroosterhoff/Documents/Werk/Brugada/DATA/Measurements/AJM060/BSPM/ajm060B.bdf';
[pathstr, name, ext] = fileparts(filename) ;
d=dir(fullfile(pathstr,[name,'+*',ext]));
fnames={d.name};
fnames=[{[name ext]},fnames]; % make sure basefile is the first

fs=2048;
flsamp=395*fs;
flsec=flsamp/fs; % file length in seconds
startsec=1060; % start ajmaline in seconds
expint=60;% interval in seconds between exports
numint=11; % number of exports including baseline
% 7 niet gedaan
for i=8:numint
    exptime=startsec+(i-1)*expint;
    fileind=1+floor(exptime/flsec); %index to file
    offset=mod(exptime,flsec); %offset within file
    prefstart=offset-4;
    prefend=offset+5;
    rstart=max(0,prefstart);
    rend=min(flsec-1,prefend);
    if prefstart<0 || prefend>flsec-1
        display(sprintf('Showing %d to %d, target at %d',rstart,rend,offset-rstart));
        pause(3);
    end
    clipboard('copy',sprintf('%ds_Ajm',(i-1)*expint));
    ViewandExport(fullfile(pathstr,fnames{fileind}),'amsbdf','1:68',num2str(rstart),num2str(rend));
end

display(['Finished ' mfilename]);  
        


