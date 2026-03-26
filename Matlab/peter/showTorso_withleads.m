global AS
global ATRIA
global TORSO
global VENTR

% PHIinit(0)

figure(1)
clf
patch('Faces',TORSO.LLITRI,'Vertices',TORSO.LLVER,...
				 'FaceLighting','phong','BackFaceLighting','unlit','AmbientStrength',0.7,...
				 'faceColor','k','edgecolor','none','FaceAlpha',0.4);
			 material([0.5 0.5 0.5])
patch('Faces',TORSO.RLITRI,'Vertices',TORSO.RLVER,...
				 'FaceLighting','phong','BackFaceLighting','unlit','AmbientStrength',0.7,...
				 'faceColor','k','edgecolor','none','FaceAlpha',0.4);
material([0.5 0.5 0.5])
patch('Faces',TORSO.TITRI,'Vertices',TORSO.TVER,...
				 'FaceLighting','phong','BackFaceLighting','unlit','AmbientStrength',0.7,...
				 'faceColor',[200 80 50]/255,'edgecolor','none','FaceAlpha',0.35);
material([0.5 0.5 0.5])
hs=patch('Faces',ATRIA.ITRI,'Vertices',ATRIA.VER,...
				 'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
				 'FaceColor','r',...
				 'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');             
			 material metal
hs=patch('Faces',VENTR.ITRI,'Vertices',VENTR.VER,...
				 'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
				 'FaceColor','r',...
				 'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');             
             material metal
material metal
             
view(90,0);axis equal off vis3d
camlight left 
camlight right
% camlight headlight

vs=[19 26 34 41 48 54];
Vlab=['V1';'V2';'V3';'V4';'V5';'V6'];
for i=1:length(vs)
	p=TORSO.TVER(vs(i),:);
	[ti,l]=find(TORSO.TITRI==vs(i));
	b=cross(TORSO.TVER(TORSO.TITRI(ti(:),1),:)-TORSO.TVER(TORSO.TITRI(ti(:),2),:),...
		  	TORSO.TVER(TORSO.TITRI(ti(:),1),:)-TORSO.TVER(TORSO.TITRI(ti(:),3),:));
	c=mean(b); c=c./norm(c);
	
	if i==3
		p=mean(TORSO.TVER(vs(i)-1:vs(i),:));
	end
	p=p+c*0.01;
	line(p(1),p(2),p(3),'Marker','o','Markerfacecolor','w','markersize',4,'color','k')
	text(p(1),p(2),p(3),Vlab(i,:),'color','k')
end

saveas(gcf,[dir 'torso' ext]);




figure(2)
clf
patch('Faces',TORSO.LLITRI,'Vertices',TORSO.LLVER,...
				 'FaceLighting','phong','BackFaceLighting','unlit','AmbientStrength',0.7,...
				 'faceColor','k','edgecolor','none','FaceAlpha',0.4);
			 material([0.5 0.5 0.5])
patch('Faces',TORSO.RLITRI,'Vertices',TORSO.RLVER,...
				 'FaceLighting','phong','BackFaceLighting','unlit','AmbientStrength',0.7,...
				 'faceColor','k','edgecolor','none','FaceAlpha',0.4);
material([0.5 0.5 0.5])
patch('Faces',TORSO.TITRI,'Vertices',TORSO.TVER,...
				 'FaceLighting','phong','BackFaceLighting','unlit','AmbientStrength',0.7,...
				 'faceColor',[200 80 50]/255,'edgecolor','none','FaceAlpha',0.35);
material([0.5 0.5 0.5])
hs=patch('Faces',ATRIA.ITRI,'Vertices',ATRIA.VER,...
				 'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
				 'FaceColor','r',...
				 'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');             
			 material metal
hs=patch('Faces',VENTR.ITRI,'Vertices',VENTR.VER,...
				 'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
				 'FaceColor','r',...
				 'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');             
             material metal
material metal
             
view(-90,0);axis equal off vis3d
camlight left 
camlight right
% camlight headlight

% vs=[19 26 34 41 48 54];
% Vlab=['V1';'V2';'V3';'V4';'V5';'V6'];
% for i=1:length(vs)
% 	p=TORSO.TVER(vs(i),:);
% 	[ti,l]=find(TORSO.TITRI==vs(i));
% 	b=cross(TORSO.TVER(TORSO.TITRI(ti(:),1),:)-TORSO.TVER(TORSO.TITRI(ti(:),2),:),...
% 		  	TORSO.TVER(TORSO.TITRI(ti(:),1),:)-TORSO.TVER(TORSO.TITRI(ti(:),3),:));
% 	c=mean(b); c=c./norm(c);
% 	
% 	if i==3
% 		p=mean(TORSO.TVER(vs(i)-1:vs(i),:));
% 	end
% 	p=p+c*0.01;
% 	line(p(1),p(2),p(3),'Marker','o','Markerfacecolor','w','markersize',4,'color','k')
% 	text(p(1),p(2),p(3),Vlab(i,:),'color','k')
% end

saveas(gcf,[dir 'torso_back' ext]);






	
