% loopnode.m
% 020226
% function loop=loopnode(VER,ITRI,node)
% finds a the loop formed by the direct neighbours of node on 
% a triangulated surface; if the loop does not include node: the loop
% indeed loops around node: size(loop)=number of triangles on which node is
% situated
% the loop is clockwise (viewed from outside) around node, 
% the same (uniform) sense is assumed for the triangulation

function loop=loopnode(VER,ITRI,node)
list=[];
[list BTRI]=buren(VER,ITRI,node);
   
if isempty(list)==0,
   [nbtri ndum]=size(BTRI);
   % cycle entries such that BTRI(:,1)=node
   for i=1:nbtri,
       trits=BTRI(i,2:4);
       l=find(trits==node);
       if l==3, trits=[trits(3) trits(1) trits(2)]; end
       if l==2, trits=[trits(2) trits(3) trits(1)]; end
       BTRI(i,2:4)=trits(1:3);
   end
  
   unpaired=[]; begin=1;
   for i=1:nbtri
      a=BTRI(i,3);
      unpaired=find(BTRI(:,4)==a);
      if isempty(unpaired)~=0, unpaired=a;, begin=i;, break, end
   end
   
   % form clockwise (viewed from outside) loop around node
   loop=BTRI(begin,2:4);
   nloop=3;
   
   nbtri=nbtri-1;
   while nbtri >0, 
      i=find(BTRI(:,3)==loop(nloop));
      if isempty(i)~=1,
         nloop=nloop+1;
         loop(nloop)=BTRI(i,4);
       end
       BTRI(i,:)=[];
       nbtri=nbtri-1;
      end
      
end
if loop(2)==loop(nloop), loop=loop(2:nloop-1); end
 
      
   
   