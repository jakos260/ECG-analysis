% refine1.m
% function [VER,ITRI]=refine1(VER,ITRI,itri)
% refine triangle ``itri'' by halving its longest edge; as well as
% the triangle that shares this edge (if present)
% if itri unspecified: treats the longest edge among ITRI
% 20050416

function [VER,ITRI]=refine1(VER,ITRI,itri)
ntri=size(ITRI,1);
nver=size(VER,1);

if nargin<3,
    % identify longest edge
     
     long=0; ilong=0;
     for i=1:ntri,
         d=norm3d(VER(ITRI(i,1),:)-VER(ITRI(i,2),:));
         if d>long, long=d; ilong=i; end
         d=norm3d(VER(ITRI(i,2),:)-VER(ITRI(i,3),:));
         if d>long, long=d; ilong=i; end
     end
     itri=ilong;
 end
 
 %line(VER(ITRI(itri,:),1),VER(ITRI(itri,:),2),VER(ITRI(itri,:),3),'col','w')

% find longest edge of triangle itri

EDGES=[];
iver=[ITRI(itri,3) ITRI(itri,1:3) ITRI(itri,1:2)];

EDGES=VER(iver(3:5),:)-VER(iver(2:4),:);
length=sum(EDGES.^2,2);
[long llong]=max(length);
edge=[iver(llong+2) iver(llong+1)]; % note: sense is reversed

% EDGES=ones(ntri,1)*edge;
% identify neighbouring triangle nbtri (sharing edge) 
btri=find(sum(ismember(ITRI,edge),2)==2);
btri(btri==itri)=[];

% refine triangle itri
nver=nver+1;
VER(nver,:)=(VER(edge(1),:)+VER(edge(2),:))/2;
ITRI(itri,:)=[iver(llong) iver(llong+1) nver];
ntri=ntri+1;
ITRI(ntri,:)=[nver iver(llong+2) iver(llong)];

if isempty(btri)==0,
   % refine neighbour triangle: btri
   iver=ITRI(btri,1:3);
   top=iver(iver~=edge(1)&iver~=edge(2));
   ITRI(btri,1:3)=[top  edge(1) nver];
   ntri=ntri+1;
   ITRI(ntri,1:3)=[nver edge(2) top];
end

