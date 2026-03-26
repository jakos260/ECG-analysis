% tricheck.m
% function [ISOLATED, ABUNDAND, duplicate,loners,sausages] = tricheck(VER,ITRI)
% 20050420
% check integrity of the triangle specs of a single triangulated closed surface
% see also trisect

function [ISOLATED, ABUNDAND, DUPLICATES,loners,sausages]=tricheck(VER,ITRI)
[nver idum]=size(VER);
[ntri jdum]=size(ITRI);

EDGES=[ITRI(:,1) ITRI(:,2);ITRI(:,2) ITRI(:,3);ITRI(:,3) ITRI(:,1)];
ind=[1:ntri 1:ntri 1:ntri];
EDGES=[EDGES ind'];
[nedge jdum]=size(EDGES);
EDGES=sortrows(EDGES,1);
EDGES=[EDGES ones(3*ntri,1)];
'      nver    ntri     nedge'
[nver ntri nedge]


'identify use of edges'
    
for i=1:nedge-1,
   if EDGES(i,4)==1,
   shared=[];   
   for j=i:nedge,
      if (EDGES(i,1)==EDGES(j,2)&EDGES(i,2)==EDGES(j,1))|...
         (EDGES(i,1)==EDGES(j,1)&EDGES(i,2)==EDGES(j,2)),
         shared=[shared j];
      end
   end
   EDGES(shared,4)=length(shared);
 end  
end   

'isolated edges:'
ISOLATED=EDGES(EDGES(:,4)==1,:); 
[nis,jdum]=size(ISOLATED);
['total: ' num2str(nis)]
if nis>0, 'node1  node2  itri', ISOLATED, end


'edges used more than twice:'
ABUNDAND=EDGES(EDGES(:,4)>2,:);
[nab,jdum]=size(ABUNDAND);
['total: ' num2str(nab)]
if nab>0,
'       node         node    triangle      number of times the edge is used'  
ABUNDAND

end



'search for duplicate triangles'
DUPLICATES=[];
for i=1:ntri-1,
indi=sort(ITRI(i,:));
   for j=i+1:ntri,
     indj=sort(ITRI(j,:));
     if indi==indj, DUPLICATES=[DUPLICATES;[i j]]; end
   end
end 
'duplicate triangles:'
DUPLICATES

next=1;
if next==1,
    'search for isolated vertices' 
  nodes=ITRI(:);
  vert=1:nver;
  iso=histc(nodes,vert);
  ISO=[iso vert'];
  'isolated vertices:'
  loners=ISO(ISO(:,1)==0,2)
end


sausages=[];
next=0;
if next==1,
  'search for sausage nodes'
   
   nsaus=0;

   for nodei=1:nver,
     buurtri=[];
     %find triangles containing nodei
     buurtri=find(ITRI(:,1)==nodei|ITRI(:,2)==nodei|ITRI(:,3)==nodei);
     % test if a circular path exists connecting all neighbour vertices
     % while not using nodei; 
     % if this is not the case: nodei is a sausage point
     nb=length(buurtri);nbb=nb;
     SELECT=ITRI(buurtri,:);
     if isempty(SELECT)==0,
        start=SELECT(1,:);
        kop=start(1); if kop==nodei, kop=start(2); end
        kop;
        SELECT(1,:)=[];
        nbb=nb-1;
        for j=2:nb,
            for k=1:nbb,
	            select=SELECT(k,:);
	            if isempty(select(find(select==kop)))==0,
	                kop=select(find(select~=kop&select~=nodei));
	                SELECT(k,:)=[];
	                nbb=nbb-1;
	                break;
                end
            end  
        end
     end
     if nbb>0,sausages=[sausages nodei]; nsaus=nsaus+1;, end 
   end
   [' found  ' num2str(nsaus) ' sausage nodes']
   if nsaus>0,sausages, end
end


