function TRIstatistics(VER,ITRI,nr)


figure (nr);clf
wd=zeros(length(VER),1);
for i=1:length(wd)
	[ti,l]=find(ITRI==i);
	b=cross(VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),2),:),VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),3),:));
	c=mean(b); c=c./norm(c); % mm
	v1=VER(i,:);
	v2=v1+c;
	TR=linetris(VER,ITRI,v1,v2);
	D=min(TR(TR(:,5)>0.1,5));
	if ~isempty(D)
		wd(i)=D;
	end
end
plot(sort(wd))
maxwd=mean(wd)+1*std(wd);
axis([1 length(wd) 0 maxwd]);
legend(['wall thickness mean ' num2str(mean(wd))]);


figure(nr+1);clf
patch ('Faces',ITRI,'Vertices',VER,'FaceVertexCData',wd,'FaceColor','interp',...
			'facealpha',1,'FaceLighting','phong','edgecolor','k','ButtonDown','selectnode');
% caxis([0.5 maxwd]);
colorbar
axis off equal vis3d

vol=calc_volume(VER,ITRI)
adj=graphdist(VER,ITRI,1);
adj(adj==0)=1000;
if (min(min(adj))<0.1)
	disp(['Edge distances smaller than 0.1 mm found, Minimum found ' num2str(min(min(adj)))]);
end
figure(nr+2);
clf;
plotAuto(VER,ITRI);


function vol=calc_volume(VER,ITRI)

vol=0;
for pp=1:length(ITRI),
	vol=vol+det3d(VER(ITRI(pp,1),:),VER(ITRI(pp,2),:),VER(ITRI(pp,3),:));
end
vol =-vol/6000;


%--------------------------------------------------------------

function plotAuto(AVER,AITRI)

disp('plot auto solid angle')
B=zeros(length(AVER));
hw=waitbar(0,'Calculating auto solid angle');
for i=1:length(AVER)
	[B(i,:),jsing]=rowforw(AVER,AITRI,AVER(i,:));
	waitbar(i/length(AVER),hw);
end
close(hw);
clf
colA=diag(B);
patch ('Faces',AITRI,'Vertices',AVER,'FaceVertexCData',colA,'FaceColor','interp',...
			'facealpha',1,'FaceLighting','phong','edgecolor','k','ButtonDown','selectnode');
axis off normal;
view(90,0); colormap(loadasci('tims.mcm'));
%caxis([-max(abs(colA)) max(abs(colA))]);
colorbar;
axis equal off vis3d
