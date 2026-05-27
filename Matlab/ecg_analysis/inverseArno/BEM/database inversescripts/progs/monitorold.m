% file monitor
% monitor results of invudl/invedls ; specify: montims
% 041023 A. van Oosterom

CMAP=loadasci('tims.mcm');
colormap(CMAP);

dim=size(ITRI);
ntri=dim(1);
xyzmin=min(VER);
xyzmax=max(VER);
len=max([xyzmax(1)-xyzmin(1) xyzmax(2)-xyzmin(2) xyzmax(3)-xyzmin(3)]);

xyzmax(:)=xyzmin(:)+len;
assen=[xyzmin(1) xyzmax(1) xyzmin(2) xyzmax(2) xyzmin(3) xyzmax(3)];
assen=assen;
assen=assen;

nrim=size(RIM,1);

rindex=[3 2 2 1 1 3];
rinzin=[1 -1 1 -1 1 -1];
for ispl=1:6,
  subplot(2,3,ispl)
  range=RANGES(rindex(ispl),:);
  k=0;
  for i=range(1):range(2),
     if zin(i)*rinzin(ispl) <= 0,
        k=k+1;
        TRI(k,:)=ITRI(i,:);
     end
  end
  if k>0,
     set(trisurf(TRI(1:k,:),VER(:,1),VER(:,2),VER(:,3),montims),'FaceColor','interp','EdgeColor','none')
  end
  view(90,0)
  axis(assen);
  axis square
  grid off
  axis off
  hold on
  if nrim > 0,
     for i=1:nrim,
       plot3([VER(RIM(i,1),1) VER(RIM(i,2),1)],[VER(RIM(i,1),2) VER(RIM(i,2),2)],...
       [VER(RIM(i,1),3) VER(RIM(i,2),3)],'-k')
     end
  end
end
hcbar=colorbar;
hcbarpos=[.92 .12 .02 .6];
set(hcbar,'position',hcbarpos);

