function save2trp(filepath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Name:    save2trp
% 
% Desc:    load a list of qtriplot commands recorded and stored by
% qtriplot.m and sav to trp, tri and mat files.
% 
% Inputs:  filepath: path to .mat file containing variable qtriplotcmds. User is presented
%                   With a file dialog if filepath is missing or empty.
% 
% Outputs: none
% 
% Created: 19-Feb-2014
% 
% Author:  Peter Oosterhoff
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('filepath','var') || isempty(filepath)
    [FileName,PathName] = uigetfile('*.mat','Select file containing qtriplot commands');
    filepath=fullfile(PathName,FileName);
end
[pathname, fname, ext] = fileparts(filepath);

load(filepath,'qtriplotcmds');

trppath=fullfile(pathname,'trp');
if ~exist(trppath,'dir')
    mkdir(trppath);
end
fnameout=fullfile(trppath,[fname '.trp']);
fid=fopen(fnameout,'w');

trinum=0;
funnum=0;



for i=1:length(qtriplotcmds)
    C=qtriplotcmds{i}; % assume 1-D cell array of 1-D cell array
    if length(C)==1 && ischar(C{1})
        fprintf(fid,'%s \n',sprintf(C{:}));
    elseif (length(C)==2 || length(C)==3) && (isnumeric(C{1})) && (isnumeric(C{2}))
        triname=sprintf('%06d.tri',trinum);
        savetri(fullfile(trppath,triname),C{1},C{2});
        trinum=trinum+1;
        if length(C)==3 && ischar(C{3})
            objstr=[C{3} '='];
        else
            objstr='';
        end
        fprintf(fid,'file %s%s\n',objstr,triname);
    elseif  (length(C)==1 || length(C)==2) && isnumeric(C{1})
        funname=sprintf('%06d.fun',funnum);
        saveasci(fullfile(trppath,funname),C{1});
        funnum=funnum+1;
        if length(C)==2 && ischar(C{2})
            objstr=[C{2} '='];
        else
            objstr='';
        end
        fprintf(fid,'funfile %s%s\n',objstr,funname);
        
    end
end
fclose(fid);

