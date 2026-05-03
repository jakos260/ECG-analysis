function [tri,edgenr]=findTriEdge(ITRI,edge)

edgenr=1;
tri=find(ITRI(:,1)==edge(1) & ITRI(:,2)==edge(2));
if isempty(tri)
	tri=find(ITRI(:,2)==edge(1) & ITRI(:,3)==edge(2));
	edgenr=2;
end
if isempty(tri)
	tri=find(ITRI(:,3)==edge(1) & ITRI(:,1)==edge(2));
	edgenr=3;
end
if isempty(tri)
	edgenr=0;
end
