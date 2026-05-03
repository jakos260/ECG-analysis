function off2tri(dirname)

files=dir('*.off');
for i=1:length(files)
    [VER,ITRI]=readoff(files(i).name);
    savetri([files(i).name(1:end-3) 'tri'],VER,ITRI);
    
end