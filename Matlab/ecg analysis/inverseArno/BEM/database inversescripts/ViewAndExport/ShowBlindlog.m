load('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Biosemi/export/Blindlog_20120702T152629.mat','blindlog');
% load('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig08/Biosemi/export/Pacing/Blindlog_20121218T202211.mat','blindlog');
for i=1:length(blindlog);
    [pathstr, name, ext] = fileparts(blindlog(i).exportfilename);
    fprintf('%d: %s\n',i,name);
end