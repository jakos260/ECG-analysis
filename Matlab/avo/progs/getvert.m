% file getvert.m
% callback of triplot

% identify the node of the set VER of triangulated surface that is
% visible on screen and nearest to the mouse pointer

L=get(gca,'currentpoint');
l1=L(1,:);
l2=L(2,:);
% intersections:
INTERS=linetris(VER,ITRI,l1,l2);
% only the facing ones are relavant
SELECT=INTERS(INTERS(:,2)>0,:);
[small,is]=min(SELECT(:,5));
iss=SELECT(is,1);
veris=ITRI(iss,:);
V1=VER(veris(1),:);
V2=VER(veris(2),:);
V3=VER(veris(3),:);
intersept=l1+SELECT(is,5)*(l2-l1);
disis(1)=(intersept-V1)*(intersept-V1)';
disis(2)=(intersept-V2)*(intersept-V2)';
disis(3)=(intersept-V3)*(intersept-V3)';
[small,ii]=min(disis);
id=veris(ii);
set(ui4,'string',sprintf('%5.3f',fun(id)));
set(ui5,'string',num2str(id));
delete(ht)
ht=text(1.02*VER(id,1),1.02*VER(id,2),1.02*VER(id,3),'0','color','w');
set(ht,'HorizontalAlignment','Center')


