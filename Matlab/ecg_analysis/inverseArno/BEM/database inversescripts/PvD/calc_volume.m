function TRIstatistics(VER,ITRI)


figure (1003);clf
wd=zeros(length(AVER),1);
for i=1:length(wd)
	[ti,l]=find(AITRI==i);
	b=cross(AVER(AITRI(ti(:),1),:)-AVER(AITRI(ti(:),2),:),AVER(AITRI(ti(:),1),:)-AVER(AITRI(ti(:),3),:));
	c=mean(b); c=c./norm(c); % mm
	v1=AVER(i,:);
	v2=v1+c;
	TR=linetris(AVER,AITRI,v1,v2);
	D=min(TR(TR(:,5)>0.1,5));
	if ~isempty(D)
		wd(i)=D;
	end
end
plot(sort(wd))

legend(['wall thickness mean ' num2str(mean(wd))]);
figure (1005);clf
patch ('Faces',AITRI,'Vertices',AVER,'FaceVertexCData',wd,'FaceColor','interp',...
			'facealpha',1,'FaceLighting','phong','edgecolor','k','ButtonDown','selectnode');
caxis([0.5 5]);
colorbar
axis off equal vis3d

vol=calc_volume(VER,ITRI)



function vol=calc_volume(VER,ITRI)

vol=0;
for pp=1:length(ITRI),
	vol=vol+det3d(VER(ITRI(pp,1),:),VER(ITRI(pp,2),:),VER(ITRI(pp,3),:));
end
vol =vol/6000;

