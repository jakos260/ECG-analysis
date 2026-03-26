% triverlabels.m
% relabels the vertices and triangle indices of a .tri file
% vertex and triangle labels are started at node, to which label 1 is
% assigned. Next, labels are assigned spiralling around this node in 
% clocksize sense as seen from the outside of the surface.

% A. van Oosterom; 2005-04-06; method based on edges

%function [VER,ITRI]=triverlabels(VER,ITRI,node); 

% clear
%[VER,ITRI]=icosasub(3);
%[VER,ITRI]=loadtri('../../broddel/uprofile.tri');
%[VER,ITRI]=loadtri('../../atria/atria_avo.tri');

node=1;

nver=size(VER,1);ntri=size(ITRI,1);
EDGES=[ITRI(:,1) ITRI(:,2);ITRI(:,2) ITRI(:,3);ITRI(:,3) ITRI(:,1)];
EDGES=sort(EDGES,2);
list=max(EDGES(:,1))*EDGES(:,1)+EDGES(:,2);
EDGES=sortrows([EDGES list]);
EDGES=skipreplica(EDGES,3);

VERSPECS=(1:nver)';
    
nedges=size(EDGES,1);   
nver=size(VER,1); ntri=size(ITRI,1);
verdist=(edgedist(nver,ITRI,node))';
EXTREMES=localextremes(ITRI,verdist);
EXTREMES=EXTREMES(EXTREMES(:,3)==1,1:2);
VERSPECS=[VERSPECS verdist];
EDGES(:,3)=mean(verdist(EDGES(:,1:2)),2);
EDGES=EDGES(EDGES(:,3)-floor(EDGES(:,3))~=0.5,:);
tridist=ceil(mean(verdist(ITRI),2));
   
TRIDIST=[(1:ntri)' tridist];
VERSPECS=[(1:nver)' verdist zeros(nver,1) zeros(nver,1)];
VERSPECS(node,3:4)=[1 0];

VALS=VERSPECS(:,2);
clf
grsw=1;
triplot

kmax=max(verdist);
list=node;
last=node;
verlab=1;
for k=1:kmax,
    stringlab=1;
    USE=EDGES(EDGES(:,3)==k,1:2);
    if isempty(USE)==1, break, end
    allverk=skipreplica(USE(:));
    nverk=length(allverk);
    % find starting point for next string
    [mi istart]=min(norm3d(VER(allverk,:)-ones(nverk,1)*VER(last,:)))
    start=allverk(istart)
    loop=findloop(ITRI,last)
    last=start;
    list=[list start]
    nlist=length(list);
    verlab=verlab+1;
    VERSPECS(list(nlist),3)=verlab;
    VERSPECS(list(nlist),4)=stringlab;
    endlist=list(nlist)
    node=endlist; 
    setnode
    % find direct neighbor of start that will setup the correct orientation of the string
    nloop=length(loop)
    ib=find(loop==start)+1;
    if ib>nloop, ib=1; end
    new=loop(ib)
    list=[list new]
    
    nlist=length(list);
    verlab=verlab+1;
    VERSPECS(nlist,3)=verlab;
    VERSPECS(nlist,4)=stringlab;
    
    USE(USE(:,1)==min(node,new)&USE(:,2)==max(node,new),:)=[]
    endlist=new;
    line(VER([node new],1), VER([node new],2), VER([node new],3),'col','w') 
    node=endlist; 
    setnode
    back=0;
    
    while isempty(USE)==0,
       connect=find(USE(:,1)==endlist|USE(:,2)==endlist)
       [endlist connect']
       ncon=size(connect,1);
       for n=1:ncon,
           if USE(connect(n),1)~=endlist, new=USE(connect(n),1); else, new=USE(connect(n),2);, end
           line(VER([endlist new],1), VER([endlist new],2), VER([endlist new],3),'col','w') 
           node=new;
           setnode;
           if any(list==new)==0,
               list=[list new];
               verlab=verlab+1;
               VERSPECS(new,3)=verlab;
               VERSPECS(new,4)=stringlab;
           end
       end
       
       USE(connect,:)=[]
       list
        if k==3, list; pause; end
       if exist('new'),endlist=new; end
       nlist=size(list,2);
       if ncon==0,
           back=back+1;
           if nlist-back>0,
               endlist=list(nlist-back); 
           else,
               endlist=USE(1,1);
               if any(list==endlist)==0,
                  verlab=verlab+1;
                  list=[list endlist];
                  VERSPECS(endlist,3)=verlab;
                  VERSPECS(endlist,4)=stringlab;
               end
               stringlab=stringlab+1;
           end
       end
   end
   list
   pause
end

% fill in any nodes in EXTREMES
if isempty(EXTREMES)==0,
   nextremes=size(EXTREMES,1);
   for j=1:nextremes,
       new=EXTREMES(j,1);
       val=EXTREMES(j,2);
       bver=buren(ITRI,new);
       stringlab=VERSPECS(bver(1),4);
       [new val-1 stringlab]
       % insert new after the last entry of [val-1; stringlab]
       last=max(VERSPECS(VERSPECS(:,2)==val-1&VERSPECS(:,4)==stringlab,3))
       
       list=VERSPECS(VERSPECS(:,2)>val-1,1)
       VERSPECS(list,3)=VERSPECS(list,3)+1;
       VERSPECS(new,3)=last+1;
   end
end


VERSPECS

'ready'

next=0;
if next==1,
   % inspect the results for induvidual k values
   for k=1:kmax,
       USE=VERSPECS(VERSPECS(:,2)==k,:);
       USE=sortrows(USE,3)
       nuse=size(USE,1);
       labmax=USE(nuse,4);
       for ll=1:labmax,
           [k ll]
           string=USE(USE(:,4)==ll,1);
           string=[string; string(1)]
           clf
           triplot
           line(VER(string,1), VER(string,2), VER(string,3),'col','w') 
           node=string(1)
           setnode
           pause
       end
   end
end
   
  next=0;
  if next==1,
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
       end
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