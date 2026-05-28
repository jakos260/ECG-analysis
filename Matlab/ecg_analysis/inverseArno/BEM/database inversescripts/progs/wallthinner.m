function wallthinner

inname='output/PPD1_h40_';
outname='output/PPD1_h40_fosa2_';
[AVER,AITRI]=loadtri([inname 'atria.tri']);
[EVER,EITRI]=loadtri([inname 'endo.tri']);
[VVER,VITRI]=loadtri([inname 'fosa_ventricle.tri']);


msk=ones(1,length(AVER));
for i=1:length(AVER);
	if ~isempty(find(EVER(:,1)==AVER(i,1) & EVER(:,2)==AVER(i,2)& EVER(:,3)==AVER(i,3)))
		msk(i)=0;
	end
end
epi=1:length(AVER);
epi=epi(find(msk));
orgAVER=AVER;
for i=2.6:0.5:9.6
	for j=1:3
		AVER=thinnerwall(AVER,AITRI,epi,i*2/3,i);
	end
	[AVER,error]=CheckGeom(AVER,AITRI,orgAVER,VVER,VITRI,EVER,EITRI);%check geometry
	if (error)
		return;
	else
		VVER=UpdateVertices(VVER,AVER,orgAVER); % 
		EVER=UpdateVertices(EVER,AVER,orgAVER); %
		orgAVER=AVER;
	end
end


figure (1003)
wd=zeros(length(AVER),1);
for i=1:length(wd)
	[ti,l]=find(AITRI==i);
	b=cross(AVER(AITRI(ti(:),1),:)-AVER(AITRI(ti(:),2),:),AVER(AITRI(ti(:),1),:)-AVER(AITRI(ti(:),3),:));
	c=mean(b); c=c./norm(c); % mm
	v1=AVER(i,:);
	v2=v1+c*10;
	TR=linetrisect(AVER,AITRI,v1,v2);
	if ~isempty(TR)
		wd(i)=10*min(TR(:,5));
	end
end
plot(sort(wd))
legend(['wall thickness mean ' num2str(mean(wd))]);



savetri([outname 'atria.tri'],AVER,AITRI);
savetri([outname 'ventricle.tri'],VVER,VITRI);
savetri([outname 'endo.tri'],EVER,EITRI);
splitendo([outname 'endo']);

plotAuto(AVER,AITRI); %plot autosolid angle
plotAuto(VVER,VITRI); %plot autosolid angle









% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VER=thinnerwall(VER,ITRI,nodes,maxdist,MAXD)

c=zeros(length(nodes),3);
for i=1:length(nodes)
	[ti,l]=find(ITRI==nodes(i));
	b=cross(VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),2),:),VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),3),:));
	c(i,:)=mean(b); c(i,:)=c(i,:)./norm(c(i,:)); % mm
end


k=0;
for i=1:length(nodes)
	[ti,l]=find(ITRI==nodes(i));
	v1=VER(nodes(i),:);
	v2=v1+c(i,:)*MAXD;
	TR=linetrisect(VER,ITRI,v1,v2);
	if ~isempty(TR)
		alpha=min(TR(:,5));
		if alpha*MAXD > maxdist
			keepVert=VER(nodes(i),:);
			VER(nodes(i),:)=VER(nodes(i),:)+c(i,:)*(alpha*MAXD-maxdist);
			tj=unique(ITRI(ti,:));
			tj(find(tj==nodes(i)))=[];
			k=k+1;
			for j=1:length(tj)
				TRT=linetrisect(VER,ITRI,VER(nodes(i),:),VER(tj(j),:));
				if ~isempty(TRT)
					keepVert;
					VER(nodes(i),:)=keepVert;  %restore due to crossection
					k=k-1;
					break;
				end
			end
		end
	end
end

disp(['adapted ' num2str(k) ' nodes from ' num2str(MAXD) ' to ' num2str(maxdist) ' mm']);

%********************************************************************

function plotAuto(AVER,AITRI)

disp('plot atria')
B=zeros(length(AVER));
for i=1:length(AVER)
	[B(i,:),jsing]=rowforw(AVER,AITRI,AVER(i,:));
end
AMA_A=B;

colA=AMA_A.*eye(size(AMA_A)); colA=colA(colA~=0);
figure(11)
clf;
patch ('Faces',AITRI,'Vertices',AVER,'FaceVertexCData',colA,'FaceColor','interp',...
			'facealpha',1,'FaceLighting','phong','edgecolor','k','ButtonDown','selectnode');
axis off normal;
view(90,0); colormap(loadasci('tims.mcm'));
%caxis([-max(abs(colA)) max(abs(colA))]);
colorbar;axis equal vis3d



%------------------------------------------------------------------

function otherVER=UpdateVertices(otherVER,VER,orgVER)
% adapt also the endocadrium, e.g. the vertices must be adapted

for i=1:length(VER)
	a=find(otherVER(:,1)==orgVER(i,1) & otherVER(:,2)==orgVER(i,2) & otherVER(:,3)==orgVER(i,3));
	if ~isempty(a)
		otherVER(a,:)=VER(i,:);
	end
end
