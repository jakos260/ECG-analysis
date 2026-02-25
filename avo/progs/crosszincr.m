% crosszincr.m
% callback of crossec.m setting step: the slider increment for zlevel 
zstep=get(ie5,'string');
zincr=str2num(zstep);
%sscanf(zstep,'%f');

slmax=ceil(slmax/zincr)*zincr;
slmin=floor(slmin/zincr)*zincr;


set(sllevel,'min',slmin,'max',slmax,'Value',zlevel);
sstep=zincr/(slmax-slmin);
slstep=[sstep .1];
set(sllevel,'Sliderstep',slstep);

