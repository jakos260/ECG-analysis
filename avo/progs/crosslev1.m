% crosslev1.m
% callback of crossec.m setting the zlevel 
% accept slider value and plot crossection
% AvO; 20050212
  zlevel=get(sllevel,'Value');
  set(ie1,'string',num2str(zlevel));
  crossub
