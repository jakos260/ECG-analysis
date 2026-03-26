% trilabels.m
% relabels the vertices and triangle indices of a .tri file
% vertex and triangle labels are started at node, to which label 1 is
% assigned. Next, labels are assigned spiralling around this node in 
% clocksize sense as seen from the outside of the surface.

% A. van Oosterom; 2005-04-03

%function [VER,ITRI]=trilabls(VER,ITRI,node); 
   clear
   %[VER,ITRI]=icosasub(1);
   [VER,ITRI]=loadtri('../atria/atriavo_1.tri');
   node=218
   nver=size(VER,1); ntri=size(ITRI,1);
   verdist=(edgedist(nver,ITRI,node))';
   tridist=ceil(mean(verdist(ITRI),2));
   
   TRIDIST=[(1:ntri)' tridist];
   VERDIST=[(1:nver)' verdist];
   
   VALS=VERDIST(:,2);
   clf
   grsw=1;
   triplot
   loop=loopnode(ITRI,node);
   nloop=length(loop);
   [mi start]=min(norm3d(VER(loop,:)-ones(nloop,1)*VER(node,:)));
   if start < nloop, 
       loop=loop([start:nloop 1:start-1]);
   else,
       loop=loop([start 1:start-1]);
   end
   loop
   last=loop(1)
   vlist=[node loop];
   line(VER([loop loop(1)],1), VER([loop loop(1)],2), VER([loop loop(1)],3),'col','w') 
   
   kmax=max(verdist);
   pause
   for k=2:kmax;
       vring=VERDIST(VERDIST(:,2)==k,1)
       tring=TRIDIST(TRIDIST(:,2)==k,1);
       nring=length(vring)
       ntring=length(tring);
       [mi start]=min(norm3d(VER(vring,:)-ones(nring,1)*VER(last,:)));
       start=vring(start);
       vlist=[vlist start]
       TEST=[ITRI(tring,:) ITRI(tring,1)]
       j=1;
       startr=start;
       while last~=startr,
           ntring=size(TEST,1);
           REF=ones(ntring,1)*[last start];
           i=find(sum(TEST(:,1:2)==REF,2)==2|sum(TEST(:,2:3)==REF,2)==2|sum(TEST(:,3:4)==REF,2)==2);
           nm=length(vlist);
           while isempty(i)==1,
               nm=nm-1;
               start=vlist(nm)
               REF=ones(ntring,1)*[last start];
               i=find(sum(TEST(:,1:2)==REF,2)==2|sum(TEST(:,2:3)==REF,2)==2|sum(TEST(:,3:4)==REF,2)==2);
           end
           trits=TEST(i,:);
           TEST(i,:)=[];
           while trits(1)~=last,
                 trits=trits([2 3 1]);
           end
           trits
           if trits(3)==startr,
               line(VER([trits(2) startr],1), VER([trits(2) startr],2), VER([trits(2) startr],3),'col','w') 
               last=startr;
               break
           end
           if any(trits(3)==vlist)==0, 
              vlist=[vlist trits(3)]
              line(VER([start trits(3)],1), VER([start trits(3)],2), VER([start trits(3)],3),'col','w') 
              start=trits(3),
          else,
              last=trits(3)
              start=trits(2)
           end
           node=start;
           setnode
       %pause
       end
   end
   
  next=0;
  if next==1,
   
 % relabel vertex indices; clockwise as seen from the outside, around node 1
   LIST=sortrows([(edgedist(nver,ITRI,1))' + (atan2(VER(:,1),VER(:,2))+pi)/(2*pi) (1:nver)'],1);
   list=[LIST(:,2) (1:nver)'];
   VER=VER(list(:,1),:);
   list=sortrows(list,1);
   ntri=size(ITRI,1);
   for i=1:ntri,
       ITRI(i,1:3)=[list(ITRI(i,1),2) list(ITRI(i,2),2) list(ITRI(i,3),2) ];
   end
   dist=(edgedist(nver,ITRI,1))';
   D=[dist(ITRI(:,1)) dist(ITRI(:,2)) dist(ITRI(:,3))];
   list=ceil(mean(D,2));
   CGRAV=(VER(ITRI(:,1),:)+VER(ITRI(:,2),:)+VER(ITRI(:,3),:))/3;
   LIST=sortrows([list + (atan2(CGRAV(:,1),CGRAV(:,2))+pi)/(2*pi) (1:ntri)'],1);
   %list=[LIST(:,2) (1:nver)'];
   %sort(LIST)
   ITRI=ITRI(LIST(:,2),:);
   
   %LIST=sortrows([(1:ntri)' sum(ITRI,2)],2);
   %ITRI=ITRI(LIST(:,1),:);
end