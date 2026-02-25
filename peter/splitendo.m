function splitendo(fname)

[VER,ITRI]=loadtri([fname '.tri']);
S=splitgraph(VER,ITRI);

VER1=VER(find(S==1),:);
VER2=VER(find(S==2),:);

ITRI1=[];
for i=1:length(S)
	if S(i)==1
		a=find(ITRI(:,1)==i);
		if ~isempty(a)
			ITRI1=[ITRI1; ITRI(a,:)];
		end
	end
end
k=0;
for i=1:length(S)
	if S(i)~=1 
		a=find(ITRI1>k);
		ITRI1(a)=ITRI1(a)-1;
	elseif S(i)==1
		k=k+1;
	end
end


ITRI2=[];
for i=1:length(S)
	if S(i)==2
		a=find(ITRI(:,1)==i);
		if ~isempty(a)
			ITRI2=[ITRI2; ITRI(a,:)];
		end
	end
end
k=0;
for i=1:length(S)
	if S(i)~=2 
		a=find(ITRI2>k);
		ITRI2(a)=ITRI2(a)-1;
	elseif S(i)==2
		k=k+1;
	end
end
clf
patch ('Faces',ITRI1,'Vertices',VER1,'FaceColor','r','facealpha',1,'FaceLighting','phong','edgecolor','none')
patch ('Faces',ITRI2,'Vertices',VER2,'FaceColor','b','facealpha',1,'FaceLighting','phong','edgecolor','none')
axis off equal
savetri([fname '_lcave.tri'],VER1,ITRI1);
savetri([fname '_rcave.tri'],VER2,ITRI2);

