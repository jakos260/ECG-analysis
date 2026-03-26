% buren.m
% function [bver,BTRI]=buren(VER,ITRI,node)
% see alsovoisin
% for node node of a triangulated
% surface: find 
%               bver=  indexes of the neighbouring vertices 
%               BTRI = triangles carrying the node

function [bver, BTRI]=buren(VER,ITRI,node)
[ntri jdum]=size(ITRI);
[nver,jdum]=size(VER);

nbver=0;
nbtri=0;
BVER=[];
BTRI=[];

index=zeros(nver,1);


% find the triangles carrying node
for i=1:ntri,
     ij=ITRI(i,:);
   if sum(ij==node)>0,
       index(ij(1:3))=ij(1:3);
       BTRI=[BTRI; i ij];
   end
end   
bver=index(index~=0);
