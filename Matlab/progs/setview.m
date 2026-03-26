% setview.m
% set the desired view on the geometry
% callback routine of triplot 
iview=get(ui7,'Value');

if ~exist('viscene'), viscene=1; end

if viscene==3,
   axes(ax1c)
else,
   axes(ax1a) 
end

%view(startpos);
if iview <=8,
   azim=90+(iview-1)*45; elev=0;
else;
   if iview==9,  azim=90; elev=90;end
   if iview==10, azim=90; elev=-90; end
end

lightpos=[cos((azim-90)/180*pi) sin((azim-90)/180*pi) 0];
set(hlight,'position',lightpos);

view(azim,elev)


