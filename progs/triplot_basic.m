% triplot_basic.m
% used to produce final figures for publications
% leaves out all uis
% NB: should be preceded by triplot; in which all parameters are set

% .....or
% specify parameters, e.g.,
% CMAP=loadmat('avojet.mcm');
% newboxsize=0.3;
% fun=zeros(size(VER,1),1);
% VALS=fun;
% contourcolor='w';
% lsw=1;
% azim=90;
% elev=0;
% symcolor=1;
% nobar=1;
% zebra=0;
% fhs=patch('Vertices',VER,'Faces',ITRI,...
% 'FaceVertexCData',fun,'FaceColor','inter');



% 20080708
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('newboxsize'),
    boxsize=max(norm3d(VER));
else,
    boxsize=newboxsize;;
end


colormap(CMAP) 
ax1a=axes('pos',[.1 .1 .85 .85]);
set(ax1a,'xlim',[-boxsize boxsize],'ylim',[ -boxsize boxsize],'zlim',[-boxsize boxsize]);
axis vis3d
axis off

% colorplot of function values involving linear interpolation
fhs=patch('Vertices',VER,'Faces',ITRI,...
'FaceVertexCData',fun,'FaceColor','inter');
set(fhs,'ButtondownFcn','nodesel');
set(fhs,'EdgeColor','none','Marker','.','MarkerEdgeColor','none');

view(azim,elev)

lightpos=[cos((azim-90)/180*pi) sin((azim-90)/180*pi) 0];
hlight=light('position',lightpos,'style','infinite');

if lsw==0,
   set(fhs ,'Facelighting','none'),
else,
   set(fhs,'facelighting','phong','Specularstrength',0.2);
end

hcbarpos=[.9 .17 .02 .7];

if symcolor==1,
   % symmetric color bar around zero
   extrmap=max(abs([fun;eps]));
   set(ax1a,'Clim',[-extrmap extrmap]);
end

if nobar==0;
   hcbar=colorbar;
   set(hcbar,'position',hcbarpos)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot contour lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('zebra')
    %if  zebra~=0,
        plotcontourlines
    %end
end

%VER=VERSAV;
% highlight selected node 
% fht=text(1.03*VER(node,1),1.03*VER(node,2),1.03*VER(node,3),'O','color','w');
% set(fht,'HorizontalAlignment','center','FontSize',16,'Fontweight','bold')

