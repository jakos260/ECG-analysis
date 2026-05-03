% addnodes.m
% function [VER,ITRI]=addnodes(VER,ITRI,NODES)
% refine triangulated surface by 
% adding a set of nodes
% refinement by selecting nearest triangle
% followed by replacing it by three subtriangles
function [VER,ITRI]=addnodes(VER,ITRI,NODES)

[nnodes,idum]=size(NODES);
for i=1:nnodes,
   pnt=NODES(i,:);
   k=isempty(find(VER(:,1)==pnt(1)&VER(:,2)==pnt(2)&VER(:,3)==pnt(3)));
   if k==1,
      % node not yet included in VER; refine
      [nver,idum]=size(VER);
      [ntri,idum]=size(ITRI);
      % find nearest triangle, based on distance to center of gravity
      large=inf;
      index=1;
      for j=1:ntri,
         ind=ITRI(j,:);
         cent=mean(VER(ind,:));
         far=norm(pnt-cent);
         if far < large, large=far; index=j; end
      end
      VER=[VER;pnt];
      nver=nver+1;
      ind=ITRI(index,:);
      [la mu dist]=lamutri(pnt,VER(ind,:),[1 2 3]);
      ITRI=[ITRI;ITRI(index,1) nver ITRI(index,3);...
                 ITRI(index,3) nver ITRI(index,2);...
                 ITRI(index,2) nver ITRI(index,1);];
      ITRI(index,:)=[];
      % test if triangulation can be improved
      
   end
end
