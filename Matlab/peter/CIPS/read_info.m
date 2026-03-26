function [INFO, info] = read_info(basedir)

basedir ='C:\Users\peter\Documents\Data\measurements';
basedir ='C:\Users\peter\Documents\Data\measurements\UCLA\prosCIPS';
info=[];
info = readDir(basedir,info);

INFO=[];
% TAB=[];
for i=1:length(info)
    if isempty(strfind(lower( info{i}.sex),'female'))
        sex=1;
    else
        sex=0;
    end
    if i==23
        a=1;
    end
    INFO=[INFO;[str2num(info{i}.id(end-1:end)) info{i}.angles info{i}.LVWallThickness info{i}.VWallVolume info{i}.circChest info{i}.age info{i}.weight info{i}.height sex]];
%     TAB(i)=[info{i}.id num2str([info{i}.angles info{i}.LVWallThickness info{i}.VWallVolume info{i}.circChest info{i}.age info{i}.weight info{i}.height]) info{i}.sex];
end


%%

function info = readDir(basedir,info)
dinfo = dir(fullfile(basedir, '*_modelinfo_extra.txt'));
if isempty(dinfo)
    ddd = dir(basedir);
    for i=1:length(ddd)
        if ~strcmp(ddd(i).name,'..') && ~strcmp(ddd(i).name,'.') && isdir(fullfile(basedir, ddd(i).name))
            info = readDir(fullfile(basedir, ddd(i).name),info);
        end
    end
    
else
    n= length(info)+1;
    info{n}.filename = basedir;
    fid=fopen(fullfile(basedir, dinfo(1).name));
    str = fgetl(fid);
    info{n}.id= str;
    str = fgetl(fid);
    info{n}.angles = str2num(str(21:end));
    str = fgetl(fid);
    info{n}.axisL= str2num(str(18:end-3));
    str = fgetl(fid);
    info{n}.LBaseCenter= str2num(str(14:end));
    str = fgetl(fid);
    info{n}.LVApex= str2num(str(8:end));
    str = fgetl(fid);
    info{n}.LVWallThickness= str2num(str(18:end));
    str = fgetl(fid);
    info{n}.VWallVolume= str2num(str(24:end));
    str = fgetl(fid);
    info{n}.AWallVolume= str2num(str(19:end));
    str = fgetl(fid);
    info{n}.circChest= str2num(str(41:end));
    str = fgetl(fid);
    info{n}.age= str2num(str(4:end-1));
    if isempty(info{n}.age)
        info{n}.age=-1;
    end
    str = fgetl(fid);
    info{n}.weight= str2num(str(7:end));
    str = fgetl(fid);
    info{n}.height= str2num(str(7:end));
    str = fgetl(fid);
    info{n}.sex= str(4:end);

    fclose(fid);
    
    
end



    
    