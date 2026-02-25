% displcol.m
% script of triplot.m 
% input and process new column
coltemp=get(sltri,'value');
column=round(coltemp);
set(sltri,'value',column);
fun=VALS(:,column);
set(hs,'FaceVertexCData',fun);


extrmap=max(abs([fun;eps]));
if symcolor==1, set(ax1a,'Clim',[-extrmap extrmap]);
   hcbar=colorbar;
   set(hcbar,'position',hcbarpos),
else,
    set(colorbar,'position',hcbarpos)
end


set(ui17,'string',num2str(column))

if exist('ui4'),
set(ui4,'string',sprintf('%3.3f',fun(node)));
end
