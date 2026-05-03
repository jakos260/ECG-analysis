function [VER,ITRI]=adaptTRI(VER,ITRI,vers);
%adaptTRI This function splits a rib between the two vertices (vers) in
%half, thus two new triangles are created. First the vertices are looked up
%if they excist in this geometry the triangle is split.


iv1=find(VER(:,1)==vers(1,1) & VER(:,2)==vers(1,2) & VER(:,3)==vers(1,3));
iv2=find(VER(:,1)==vers(2,1) & VER(:,2)==vers(2,2) & VER(:,3)==vers(2,3));
if isempty(iv1) || isempty(iv2)
	return
end
if length(iv1) > 1
	iv1=iv1(1);
end
if length(iv2) > 1
	iv2=iv2(1);
end
newv=vers(1,:)+(vers(2,:)-vers(1,:))./2; 
if ~isempty(find(VER(:,1)==newv(1,1) & VER(:,2)==newv(1,2) & VER(:,3)==newv(1,3),1))
	return
end
niv=length(VER)+1;

addtri=[];
deltri=[];
for i =1:2
	if i==1
		[tri,edgenr]=findTriEdge(ITRI,[iv1 iv2] );
	else
		[tri,edgenr]=findTriEdge(ITRI,[iv2 iv1] );
	end
	if length(tri)>1
		tri=tri(1);
	end
% 	if edgenr==0
% 		return;
% 	else
	if edgenr==1
		addtri=[addtri;...
				ITRI(tri,1) niv ITRI(tri,3);...
				ITRI(tri,3) niv ITRI(tri,2)];
		deltri=[deltri tri];
	elseif edgenr==2
		addtri=[addtri;...
				ITRI(tri,1) ITRI(tri,2) niv ;...
				ITRI(tri,1) niv ITRI(tri,3)];
		deltri=[deltri tri];
	elseif edgenr==3
		addtri=[addtri;...
			ITRI(tri,1) ITRI(tri,2) niv ;...
			ITRI(tri,3) niv ITRI(tri,2)];
		deltri=[deltri tri];
	end		
end
if isempty(addtri)
	return;
end
deltri=sort([deltri]);
for i=length(deltri):-1:1
	ITRI=[ITRI(1:deltri(i)-1,:);ITRI(deltri(i)+1:end,:)];  %remove the triangle
end
VER=[VER;newv];
ITRI=[ITRI; addtri];	


