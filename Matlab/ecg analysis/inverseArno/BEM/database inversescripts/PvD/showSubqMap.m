function showSubqMap(VER,ITRI,A,sym,minmax,delta,docolorbar)

if sym
	absmax=abs(minmax);
	minV=-absmax;
	maxV=absmax;
else
	minV=minmax(1);
	maxV=minmax(2);
end

patch('Faces',ITRI,'Vertices',VER,...
	  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
	  'FaceVertexCData',A,'FaceColor','interp',...
	  'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');
view(0,-90);
axis off equal tight
caxis([minV maxV]);
contourlines(VER,ITRI,A,'sym',1,'extremes',1,'delta',delta);             
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
if docolorbar
	b=colorbar('YTickLabel',at,'YTick',a);
	width=get(b,'TickLength');
	set(b,'TickLength',width*4);
end



