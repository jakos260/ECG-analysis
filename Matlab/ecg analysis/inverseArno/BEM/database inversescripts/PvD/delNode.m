function [VER,ITRI]=delNode(VER,ITRI,node)

VER(node,:)=[];
a=find(ITRI(:,1)~=node & ITRI(:,2)~=node &ITRI(:,3)~=node);
ITRI=ITRI(a,:);
a=find(ITRI(:,1)>node);
ITRI(a,1)=ITRI(a,1)-1;
a=find(ITRI(:,2)>node);
ITRI(a,2)=ITRI(a,2)-1;
a=find(ITRI(:,3)>node);
ITRI(a,3)=ITRI(a,3)-1;






