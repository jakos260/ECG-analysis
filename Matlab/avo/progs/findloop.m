% findloop.m
% function [loopnodes,ntrisn,nln]=findloop(ITRI,node)
% loopnodes: the direct neighbours of node on a triangulated surface
% ntris: the number of triangles of which node is a vertex
% nln:   the number of nodes of the loop; 
%        for a closed loop: nln=ntrisn and loopnodes are ordered
%        in the sense of the triangles indexes, which is assumed to be
%        identical for all triangles
%        else, nln is zero, i.e., loop is not closed
%        
% see also: loopnode, findedge, buren, voisin

% A. van Oosterom; 20050406

function [loopnodes,ntrisn,nln]=findloop(ITRI,node)
ntrisn=0; nln=0;
loopnodes=[];
[bver,BBB]=buren(ITRI,node);
list=find(ITRI(:,1)==node|ITRI(:,2)==node|ITRI(:,3)==node);
if isempty(bver)==1, return, end

% specs of triangles carrying node:
BTRI=BBB(:,2:4);
ntrisn=size(BTRI,1);
nbtri=ntrisn;

% form loop around node
trits=BTRI(nbtri,:);
l=find(trits==node);
if l==3, trits=[trits(3) trits(1) trits(2)]; end
if l==2, trits=[trits(2) trits(3) trits(1)]; end
loopnodes(1:2)=trits(2:3);
nln=2;
BTRI(nbtri,:)=[];
nbtri=nbtri-1;

while nbtri >0,
    k=[];
   [i,j]=find(BTRI==loopnodes(nln));
   if isempty(i)~=1,
        [k]=find(BTRI(i(1),:)~=loopnodes(nln)&BTRI(i(1),:)~=node);
   end     
   if isempty(k)~=1,
       nln=nln+1; 
       kk=k(1);
       [node nln kk BTRI(i(1),kk)];
       loopnodes(nln)=BTRI(i(1),kk);
       BTRI(i(1),:)=[];
       nbtri=nbtri-1;
    else,
       nln=0; break, end  
end

if nln==0, return, end

if loopnodes(nln)==loopnodes(1), 
   loopnodes(nln)=[];
   nln=nln-1;
end
 
   