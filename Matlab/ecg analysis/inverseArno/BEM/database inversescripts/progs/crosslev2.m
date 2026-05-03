% setlev2.m
% callback od crossec.m setting the level 
% accept slider value and plot crossection
% AvO; 20050212

newlev=get(ie1,'string');
zlevel=str2num(newlev);
zlevel=min(zlevel,slmax);
zlevel=max(zlevel,slmin);
set(ie1,'string',num2str(zlevel));
set(sllevel,'val',zlevel)
crossub
