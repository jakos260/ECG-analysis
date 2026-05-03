function showTorsoMap(ni,A,p,sym,minmax)

global TORSO
docont=0;
if p>0
if ni >0
	A=A-ones(length(TORSO.TVER),1)*A(ni,:);
elseif ni==0
	wct=[141 103 143];
	A=A-ones(length(TORSO.TVER),1)*mean(A(wct,:));
end
A=A(:,p);
end
clf;
if sym
	absmax=min(max(abs(A)),abs(minmax));
	A(A>absmax)=absmax;
	A(A<-absmax)=-absmax;
else
	maxV=max(max(A));maxV=max(maxV,minmax(2));
	minV=min(min(A));minV=min(minV,minmax(1));
	A(A>maxV)=maxV;
	A(A<minV)=minV;
	
end

axes('Position',[0.01 0.01 0.45 0.9])	
patch('Faces',TORSO.TITRI,'Vertices',TORSO.TVER,...
	  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
	  'FaceVertexCData',A,'FaceColor','interp',...
	  'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');
view(90,0);
axis off equal tight
if sym
	absmax=min(max(abs(A)),abs(minmax));
	caxis([-absmax absmax])
else
	caxis([minV maxV])
end
if docont
	contourlines(TORSO.TVER,TORSO.TITRI,A,'sym',sym,'extremes',1);             
end
if ni>0
	line(TORSO.TVER(ni,1),TORSO.TVER(ni,2),TORSO.TVER(ni,3),'Marker','o','color','w','markerfacecolor','w')
	distline(TORSO.TVER,TORSO.TITRI,ni,45,TORSO.DIS)
end

vs=[19 26 34 41 48 54];
Vlab=['V1';'V2';'V3';'V4';'V5';'V6'];
for i=1:length(vs)
	p=TORSO.TVER(vs(i),:);
	[ti,l]=find(TORSO.TITRI==vs(i));
	b=cross(TORSO.TVER(TORSO.TITRI(ti(:),1),:)-TORSO.TVER(TORSO.TITRI(ti(:),2),:),...
		  	TORSO.TVER(TORSO.TITRI(ti(:),1),:)-TORSO.TVER(TORSO.TITRI(ti(:),3),:));
	c=mean(b); c=c./norm(c);
		c=get(gca,'CameraTarget')-get(gca,'CameraPosition');c=-11*c/norm(c);
	if i==3
		p=mean(TORSO.TVER(vs(i)-1:vs(i),:));
		c=c*10;
	end
	p=p+c*0.05;
	line(p(1),p(2),p(3),'Marker','o','Markerfacecolor','w','markersize',4,'color','k')
% 	text(p(1),p(2),p(3),Vlab(i,:),'color','k')
end

axes('Position',[0.53 0.01 0.45 0.9])
patch('Faces',TORSO.TITRI,'Vertices',TORSO.TVER,...
	  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
	  'FaceVertexCData',A,'FaceColor','interp',...
	  'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');
view(-90,0);
axis off equal tight
if sym
	absmax=min(max(abs(A)),abs(minmax));
	caxis([-absmax absmax])
else
	minV=max(min(A),minmax(1));
	maxV=min(max(A),minmax(2));
	caxis([minV maxV])
end
if ni>0
	line(TORSO.TVER(ni,1),TORSO.TVER(ni,2),TORSO.TVER(ni,3),'Marker','o','color','w','markerfacecolor','w')
	distline(TORSO.TVER,TORSO.TITRI,ni,45,TORSO.DIS)
end
if docont
	contourlines(TORSO.TVER,TORSO.TITRI,A,'sym',sym,'extremes',1);             
end
axes('Position',[0.4 0.1 0.07 0.8])	
axis off
if sym
	absmax=min(max(abs(A)),abs(minmax));
	caxis([-absmax absmax])
else
	minV=max(min(A),minmax(1));
	maxV=min(max(A),minmax(2));
	caxis([minV maxV])
end

colorbar
		end
% if (ni)
% 	line(TORSO.TVER(ni,1),TORSO.TVER(ni,2),TORSO.TVER(ni,3),'Marker','o','color','w','markerfacecolor','w')
% end
