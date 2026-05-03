function RIM=getrim(ITRI,z,r)
% compute rim of (a part of) a triangulated surface (triangles ITRI)
% z is the sense of triangle indexes;
% r is the range of the triangles 

% get all edges

EDGES=[];
EDGES=[[min(ITRI(r(1):r(2),1:2)'); max(ITRI(r(1):r(2),1:2)'); z(r(1):r(2))']';...
       [min(ITRI(r(1):r(2),2:3)'); max(ITRI(r(1):r(2),2:3)') ;z(r(1):r(2))']';...
       [min([ITRI(r(1):r(2),3) ITRI(r(1):r(2),1)]'); max([ITRI(r(1):r(2),3)...
       ITRI(r(1):r(2),1)]');...
       ;z(r(1):r(2))']'];
EDGES=sortrows(EDGES,[1 2]);

%rim spotting 
nrim=0;
dim=size(EDGES);
n=dim(1);
i=1;
while i <= n-1,
  if EDGES(i,1:2)==EDGES(i+1,1:2),
     i=i+2;
     if EDGES(i-2,3)~=EDGES(i-1,3),
       nrim=nrim+1;
       RIM(nrim,:)=[EDGES(i-2,1) EDGES(i-2,2)];
     end
  else 
    nrim=nrim+1;
    RIM(nrim,:)=[EDGES(i,1) EDGES(i,2)];
    i=i+1;
  end
end
if i==n,
  nrim=nrim+1;
  RIM(nrim,:)=[EDGES(i,1) EDGES(i,2)];
end
