% file nodesel.m
% identify the node of the set VER of triangulated surface that is
% visible on screen and nearest to the mouse pointer
% A. van Oosterom
% bug repaired  2011_01_28:   (INTERS identified by linetris includes
% intersections outside the triangle for extremely small values of
% INTERS(:,5))
% see setnode for monitoring node signal

L=get(gca,'currentpoint');
l1=L(1,:);
l2=L(2,:);

% intersections:
INTERS=linetris(VER,ITRI,l1,l2);

INTERS(INTERS(:,3)<0 | INTERS(:,3)>1 | INTERS(:,4)<0 | INTERS(:,4)>1,:)=[];    % ###############################

% use next line if only the facing triangles are relevant
% SELECT=INTERS(INTERS(:,2)>0,:);

SELECT=INTERS;

[small,is]=min(SELECT(:,5));
iss=INTERS(is,1);

veris=ITRI(iss,:);
v1=VER(veris(1),:);
v2=VER(veris(2),:);
v3=VER(veris(3),:);
intersept=l1+SELECT(is,5)*(l2-l1);
disis(1)=(intersept-v1)*(intersept-v1)';
disis(2)=(intersept-v2)*(intersept-v2)';
disis(3)=(intersept-v3)*(intersept-v3)';
[small,ii]=min(disis);
node=veris(ii);
setnode
