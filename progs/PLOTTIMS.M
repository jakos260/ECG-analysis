function plottims(tims,geom)
[VER,ITRI]=loadtri(geom);
CMAP=loadasci('tims.mcmap');
CMAP=CMAP/255;
subplot(1,2,1)
set(trisurf(ITRI,VER(:,1),VER(:,2),VER(:,3),tims),'FaceColor','interp','EdgeColor',[.6 .6 .6])
xlabel('x-axis');
ylabel('y-axis');
colormap(CMAP)
view(90,0)
axis square
axis([-.05 .15  -.035 .105 -.07 .06])
subplot(1,2,2)
set(trisurf(ITRI,VER(:,1),VER(:,2),VER(:,3),tims),'FaceColor','interp','EdgeColor',[.6 .6 .6])
xlabel('x-axis');
ylabel('y-axis');
colormap(CMAP)
view(-90,0)
axis square
axis([-.05 .15  -.035 .105 -.07 .06])

