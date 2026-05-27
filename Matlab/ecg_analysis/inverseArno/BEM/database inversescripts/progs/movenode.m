% movenode.
% date: 20090120

figure(2)
[x,y,button]=ginput(2);
button

next=1;
if next==1,
    % process nodes of VER
'select position of a node and then its new position'

nver=size(VER,1);

V=ones(nver,1)*[-y(1) x(1) zlevel];
d=norm3d(VER-V);
[mi imi]=min(d);


VER(imi,1:2)=[-y(2) x(2) ];
VERA(imi,1:2)=[-y(2) x(2) ];

end



next=0;
if next==1,
    % process nodes of SELPNTS
'select position of a node and then its new position'
% find nearest node of VER and move it to the second position
nvers=size(SELPNTS,1);

V=ones(nvers,1)*[-y(1) x(1) zlevel];
d=norm3d(SELPNTS-V);
[mi imi]=min(d);
SELPNTS(imi,1:2)=[-y(2) x(2)];
VER(track(imi),1:2)=SELPNTS(imi,1:2);
end

figure(1)
clf
triplot
% hold on
% line(VER(track,1),VER(track,2),VER(track,3),'color','w','linewidth',1.5)
% text(VER(track(imi),1),VER(track(imi),2),VER(track(imi),3),'*','color','k','fontsize',16)  

figure(2)
ie6val=1;
crossec

