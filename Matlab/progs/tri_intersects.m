% tri_intersects.m
% function SECTS=tri_intersects(VER1,ITRI1,VER2,ITRI2);
% lists, if present, the intersections of the edges of all triangles
% of a triangulated surface S1 with any of the triangles of S2
% SECTS(:,1:6)=[ edgenode(1) edgenode(2) (on S1) triangle index of S2 of
% itersection and its verices
% NOTE 1: one way only!; 
% for autosect: just include duplicate geometries
% intersection at the node

% 20090330; A. van Oosterom

function SECTS=tri_intersects(VER1,ITRI1,VER2,ITRI2)

auto=0;
  nver1=size(VER1,1);
  ntri1=size(ITRI1,1);
  nver2=size(VER2,1);
  ntri2=size(ITRI2,1);
  if nver1==nver2,
     if sum(sum(abs(VER1-VER2)))==0,
         auto=1;
     end
  end
  
  EDGES=[ITRI1(:,1) ITRI1(:,2);
         ITRI1(:,2) ITRI1(:,3);
         ITRI1(:,3) ITRI1(:,1)];
 LIST=sort(EDGES,2);
 
SINGLES=[];
while isempty(LIST)==0,
    SINGLES=[SINGLES;LIST(1,:)];
    k=find(sum(LIST==ones(size(LIST,1),1)*LIST(1,:),2)==2);
    if isempty(k)==0,
        LIST(k,:)=[];
    end
end

EDGES=SINGLES;
  
  nedges=size(EDGES,1);
  %del=eps;
  del=sqrt(eps);
  del=del*max(max(VER1)-min(VER1));
  
  
  SECTS=[]; 
  for i=1:nedges,
      TRISECS=linetris(VER2,ITRI2,VER1(EDGES(i,1),:),VER1(EDGES(i,2),:));
      if isempty(TRISECS)==0,
         nn=size(TRISECS,1);
         SELECT=[(1:nn)' (TRISECS(:,3:5) > del) (TRISECS(:,3:5) < 1-del)];
         select=SELECT(all(SELECT(:,2:7)==1,2),1);
         nsel=length(select);
         SECTS=[SECTS; ones(nsel,1)*EDGES(i,:) TRISECS(select,1) ITRI2(TRISECS(select,1),:)];    
      end              
  end
  

  if isempty(SECTS)==0,
     if auto==1,
     % discard any (auto)intersecions at the vertices of the reported
     % triangles
      SECTS(sum(ismember(SECTS(:,[1 2]),ITRI1(SECTS(:,3),:)),2)>0    ,:)=[];
     end
  end
   
  
  
if isempty(SECTS)==0,
    'edge from  node to node    crosses TRI ; having vertices:' 
    '         of geom 1                    of geom 2  '   
end

    
