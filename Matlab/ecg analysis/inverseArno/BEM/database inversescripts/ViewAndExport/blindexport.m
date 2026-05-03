% blindexport.m
% Copy selected .mat export files to subfolder blind
% Remove comments and filename that could identify the recording
% Removed information is saved in Blindlog*.mat
clear all
clear global
close all
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select STW export files.'...
    ,'/Users/peteroosterhoff/Documents/Werk/STW/Data/','MultiSelect','on');
L=length(FileName);
if L>0
    [d,IX]=sort(rand(L,1));
    mkdir(PathName,'blind');
    dstr=datestr(now,30);
    for i=1:L
        matpath=fullfile(PathName,FileName{IX(i)});
        clear D DATA
        load(matpath,'D','DATA');
        blindlog(i).exportfilename=matpath;
        blindlog(i).filepath=D.filepath;
        D.filepath='';
        blindlog(i).comment=D.comment;
%         D.comment='*********';
        blindpath=fullfile(PathName, 'blind' ...
            ,['Blind_'  dstr '_File' num2str(i,'%02d')]);
        save(blindpath,'D','DATA');
        
    end
    save(fullfile(PathName,['Blindlog_' dstr]),'blindlog');
end