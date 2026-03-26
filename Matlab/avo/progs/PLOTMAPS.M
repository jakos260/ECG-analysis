function plotmaps(phi,geom,tekst)
[VER,ITRI]=loadtri(geom);
dim=size(VER);
nver=dim(1);
CMAP=loadasci('avo.mcmap');
CMAP=CMAP/255;

data=phi;

subplot('position',[0.0 0.2 .50 .50])
set(trisurf(ITRI,VER(:,1),VER(:,2),VER(:,3),data),'FaceColor','interp','EdgeColor','none');
colormap(CMAP);
view(90,0);
axis off
axis normal
text(.0, -.1 ,.25 ,tekst)

subplot('position',[0.5 0.2 .5 .5])
set(trisurf(ITRI,VER(:,1),VER(:,2),VER(:,3),data),'FaceColor','interp','EdgeColor','none');
colormap(CMAP);
view(270,0);
%axis equal
%axis vis3d
axis off
axis normal


