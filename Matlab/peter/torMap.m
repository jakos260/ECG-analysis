function torMap(ni,PHI,av,depV,dir ,fn)

global TORSO
global AS

if ni >0
	PHI_r=PHI-ones(length(TORSO.TVER),1)*PHI(ni,:);
elseif ni==0
	wct=[141 103 143];
	PHI_r=PHI-ones(length(TORSO.TVER),1)*mean(PHI(wct,:));
end
% PHI_r=PHI;
infor1=GetAmplInfos(PHI_r,max(depV)+av );

hf=figure(2031);set(hf,'PaperPositionMode','auto')
set(gcf,'position',[520   100   1011 406  ]);clf; 
colormap(AS.AVO)

CP=infor1(:,1);
CQ=infor1(:,2);
CT=infor1(:,3);
CST=infor1(:,4);
axes('Position',[1-1/3*(3)+0.0 0.01 0.3 0.9 ]); showMap(ni,CP,1,0.4,0.1)
axes('Position',[1-1/3*(2)+0.0 0.01 0.3 0.9 ]); showMap(ni,CQ,1,4,1)
axes('Position',[1-1/3*(1)+0.0 0.01 0.3 0.9 ]); showMap(ni,CT,1,0.8,0.25)
if ni >0
	saveas(gcf,[dir 'SurfAmplpots ' num2str(ni) ' ' fn]); 
elseif ni==0
	saveas(gcf,[dir 'SurfAmplpots_WCT ' fn]); 	
else
	saveas(gcf,[dir 'SurfAmplpots_mean ' fn]); 	
end

hf=figure(2032);set(hf,'PaperPositionMode','auto')
QRSint=sum((PHI_r(:,av:floor(av+max(depV)))),2);
set(gcf,'position',[520   100   1011 406  ]);clf; 
colormap(AS.AVO);%(204:end,:))
showMapST(ni,QRSint,1,max(abs(QRSint)),25)
if ni >0
	saveas(gcf,[dir 'QRSint ' num2str(ni) ' ' fn]); 
elseif ni==0
	saveas(gcf,[dir 'QRSint_WCT ' fn]); 	
else
	saveas(gcf,[dir 'QRSint_mean ' fn]); 
end


hf=figure(2033);set(hf,'PaperPositionMode','auto')
set(gcf,'position',[520   100   1011 406  ]);clf; 
colormap(AS.AVO)
showMapST(ni,CST,1,5,0.25)
if ni >0
	saveas(gcf,[dir 'SurfST ' num2str(ni) ' ' fn]); 
elseif ni==0
	saveas(gcf,[dir 'SurfST_WCT ' fn]); 	
else
	saveas(gcf,[dir 'SurfST_mean ' fn]); 
end



%%-----------------------------------------------------------

function showMap(ni,A,sym,minmax,delta)

global TORSO
global AS


if sym
	absmax=min(max(abs(A)),abs(minmax));
	A(A>absmax)=absmax;
	A(A<-absmax)=-absmax;
else
	A(A>maxV)=maxV;
	A(A<minV)=minV;
	
end

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
	minV=max(min(A),minmax(1));
	maxV=min(max(A),minmax(2));
	caxis([minV maxV])
end
contourlines(TORSO.TVER,TORSO.TITRI,A,'sym',sym,'extremes',1,'delta',delta);             

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
	p=p+c*0.03;
	line(p(1),p(2),p(3),'Marker','o','Markerfacecolor','w','markersize',4,'color','k')
% 	text(p(1),p(2),p(3),Vlab(i,:),'color','k')
end
colorbar
if ni>0
	distline(TORSO.TVER,TORSO.TITRI,ni,45,TORSO.DIS)
end


%%-----------------------------------------------------------
function showMapST(ni,A,sym,minmax,delta)

global TORSO
global AS


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

axes('Position',[0.01 0.01 0.48 0.9])	
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
contourlines(TORSO.TVER,TORSO.TITRI,A,'sym',sym,'extremes',1,'delta',delta);             
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

axes('Position',[0.51 0.01 0.48 0.9])	
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
% 	line(TORSO.TVER(ni,1),TORSO.TVER(ni,2),TORSO.TVER(ni,3),'Marker','o','color','w','markerfacecolor','w')
	distline(TORSO.TVER,TORSO.TITRI,ni,45,TORSO.DIS)
end
contourlines(TORSO.TVER,TORSO.TITRI,A,'sym',sym,'extremes',1,'delta',delta);             

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
% if (ni)
% 	line(TORSO.TVER(ni,1),TORSO.TVER(ni,2),TORSO.TVER(ni,3),'Marker','o','color','w','markerfacecolor','w')
% end
