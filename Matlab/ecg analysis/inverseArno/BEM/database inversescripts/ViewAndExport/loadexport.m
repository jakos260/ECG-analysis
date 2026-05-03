% loadexport: Read STW .mat export file and place begin marker at selected
% beat
function loadexport(varagin)
global bs_rejected bs_Begins
filepath='';
if nargin==1
    filepath=varagin;
    [pathstr, name, ext] = fileparts(filepath)
    [PathName,name,ext] = fileparts(filepath);
    Filename=[name '.' ext];
    FilterIndex=1;
end
if isempty(filepath)
    [FileName,PathName,FilterIndex] = uigetfile('*.mat','Select STW export file.','/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig03/BasketEnSok/export/');
    filepath=fullfile(PathName,FileName);
    [pathstr, name, ext] = fileparts(filepath)
end
clipboard('copy',name);
oldrejected=bs_rejected;
load(filepath);
display(filepath);
if ~isfield(D,'filepath')
    D.filepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig03/BasketEnSok/varken2443_3BasketEnSok+002.bdf';
end
bs_load(D.filepath,'FirstSample',D.startpos*2048+1,'lastsample',(D.endpos+1)*2048);
bs_Begins=DATA.BEATS(find(DATA.SELBEATS));
bs_rejected=oldrejected;