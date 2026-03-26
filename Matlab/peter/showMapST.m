function showMapST(ni,A,sym,minmax,delta)

global TORSO

if sym
	absmax=abs(minmax);
	minV=-absmax;
	maxV=absmax;
else
	minV=minmax(1);
	maxV=minmax(2);
end

axes('Position',[0.01 0.01 0.48 0.9])	
patch('Faces',TORSO.TITRI,'Vertices',TORSO.TVER,...
	  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
	  'FaceVertexCData',A,'FaceColor','interp',...
	  'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');
view(90,0);
axis off equal tight
caxis([minV maxV])
contourlines(TORSO.TVER,TORSO.TITRI,A,'sym',1,'extremes',1,'delta',delta);             

vs=[19 26 34 41 48 54];
% Vlab=['V1';'V2';'V3';'V4';'V5';'V6'];
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
	p=p+c*0.03;
	line(p(1),p(2),p(3),'Marker','o','Markerfacecolor','w','markersize',4,'color','k')
% 	text(p(1),p(2),p(3),Vlab(i,:),'color','k')
end
if ni>0
	line(TORSO.TVER(ni,1),TORSO.TVER(ni,2),TORSO.TVER(ni,3),'Marker','o','color','w','markerfacecolor','w')
	distline(TORSO.TVER,TORSO.TITRI,ni,45,TORSO.DIS)
end


axes('Position',[0.51 0.01 0.48 0.9])	
patch('Faces',TORSO.TITRI,'Vertices',TORSO.TVER,...
	  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
	  'FaceVertexCData',A,'FaceColor','interp',...
	  'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');
view(-90,0);
axis off equal tight
caxis([minV maxV])
contourlines(TORSO.TVER,TORSO.TITRI,A,'sym',1,'extremes',1,'delta',delta);             

width=0.07;
axes('Position',[0.4 0.1 width 0.8])	
axis off
caxis([minV maxV]);


if minV<0
	a1=[-delta:-delta:minV];a1=a1(end:-1:1)';at1=num2str(a1,'%+6.3g');
	spac='';for i=1:size(at1,2), spac=[' ' spac];  end
	if size(at1,1) >=5,for i=size(at1,1):-2:1,	at1(i,:)=spac;end;end
else
	a1=[];at1=[];
end
if maxV>0
	a2=[delta:delta:maxV]';at2=num2str(a2,'%+6.3g');
	spac='';for i=1:size(at2,2), spac=[' ' spac];  end
	if size(at2,1) >=5,for i=1:2:size(at2,1),at2(i,:)=spac;end;end
else
	a2=[];
	at2=[];
end
at0=num2str(0,'% 6.3g');for i=length(at0)+1:max(size(at1,2),size(at2,2)), at0=[' ' at0 ];  end
at=[at1; at0; at2];
a=[a1; 0; a2];
at=[at2; at0; at1];a=[a2; 0; a1];
colorbar('TickLength',[width/4 width/4],'YTickLabel',at,'YTick',a);
% line(TORSO.TVER(ni,1),TORSO.TVER(ni,2),TORSO.TVER(ni,3),'Marker','o','color','w','markerfacecolor','w')
