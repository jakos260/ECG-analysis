% autosect.m
% 20040803
% function autosect(VER,ITRI);
% lists, if any, all intersections between the edges of all triangles
% of a triangulated surface and any of the non-neighbouring
% triangles
% SECTS(:,1:3)=[ edgenode(1) edgenode(2) triangle index of itersection

function SECTS=autosect(VER,ITRI)
  [nver idum]=size(VER);
  [ntri jdum]=size(ITRI);
  SECTS=[];
  EDGES=[ITRI(:,1) ITRI(:,2);ITRI(:,2) ITRI(:,3);ITRI(:,3) ITRI(:,1)];
  EDGES=[min(EDGES'); max(EDGES')]';
  EDGES=sortrows(EDGES);
  [nedges jdum]=size(EDGES);
  list=1;
  for i=2:nedges,
      if (EDGES(i,1)~= EDGES(i-1,1))|(EDGES(i,2)~=EDGES(i-1,2)), list=[list i]; end
  end
  EDGES=EDGES(list,:);
  [nedges jdum]=size(EDGES);

  del=sqrt(eps);
  for i=1:nedges,
      TRIS=linetris(VER,ITRI,VER(EDGES(i,1),:),VER(EDGES(i,2),:));
      if isempty(TRIS)==0,
         [nn jdum]=size(TRIS);
         for j=1:nn,
             if   (TRIS(j,3) > del) & (TRIS(j,3) < 1-del)& ...
                  (TRIS(j,4) > del) & (TRIS(j,4) < 1-del)& ...
                  (TRIS(j,5) > del) & (TRIS(j,5) < 1-del),
                  SECTS=[SECTS; [EDGES(i,1:2) TRIS(j,1)]];
             end
         end
      end              
  end
if isempty(SECTS)==0,
    'edge from node: to node: crosses triangle:'
end
    
