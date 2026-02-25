function showDepRep(varargin)


node=0;onode=0;
if length(varargin)<3
	error('Not enough parameters!!!');
else
	VER=varargin{1};
	ITRI=varargin{2};
	dep=varargin{3};
	rep=varargin{4};
	maxdep=max(dep);
	mindep=min(dep);
	maxrep=max(rep);
	minrep=min(rep);
	lines=[];posi=[];
	pp=5;
	while pp<=length(varargin)
		if ischar(varargin{pp})
			key=lower(varargin{pp});
			switch key
				case 'nodes'
					node=varargin{pp+1};pp=pp+2;
				case 'onodes'
					onode=varargin{pp+1};pp=pp+2;
				otherwise
					error('unknown parameter');
			end				
		elseif size(varargin{pp},2)==3 % colormap
			lines=varargin{pp};	pp=pp+1;
		elseif length(varargin{pp})==1
			node=varargin{pp}; pp=pp+1;
			if length(varargin) > pp,
				onode=varargin{pp}; pp=pp+1;
			end
		elseif length(varargin{pp})==2
			maxdep=varargin{pp};pp=pp+1;
			mindep=maxdep(1);
			maxdep=maxdep(2);
		elseif length(varargin{pp})==4
			posi=varargin{pp};pp=pp+1;
		else
			disp(['parameter ' num2str(pp) 'not used']); 
		end			
	end
end


if isempty(lines)
	lines=loadmat('tims.mcm');
end
if abs(mindep)>10
	mindep=min(10*floor(mindep/10));
end
colormap(lines);
if isempty(posi)
	clf
	pos=get(gcf,'Position');
	set(gcf,'Position',[pos(1)  pos(2) 800   250]);
	posi=[0 0 1 1];
end

%%
axes('Position',[posi(1) posi(2) 0.55*posi(3) 1*posi(4)]);
patch('Faces',ITRI,'Vertices',VER,...
	  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
	  'FaceVertexCData',dep,'FaceColor','interp',...
	  'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');
axis off equal; 
% caxis([min(10*floor(dep/10)),maxdep]);
caxis([mindep,maxdep]);
view(90,0);
material([0.5 0.5 0.5])
camlight left
if posi(3)>.5
	if (node>0),line(VER(node,1),VER(node,2),VER(node,3),'marker','o','markersize',8,'color','w','MarkerFaceColor','w','linestyle','none'); end
	if (onode>0),line(VER(onode,1),VER(onode,2),VER(onode,3),'marker','s','markersize',8,'color','w','MarkerFaceColor','w','linestyle','none'); end
	contourlines(VER,ITRI,dep,'delta',10,'labels',0)  
else
	contourlines(VER,ITRI,dep,'delta',20,'labels',0)  
end
axes('Position',[posi(1)+posi(3)*0.35 posi(2)+posi(4)*0.2 0.1*posi(3) 0.6*posi(4)]); axis off
caxis([mindep maxdep])
colorbar;

%%
axes('Position',[posi(1)+0.49*posi(3) posi(2) 0.55*posi(3) 1*posi(4)]);
patch('Faces',ITRI,'Vertices',VER,...
	  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
	  'FaceVertexCData',rep,'FaceColor','interp',...
	  'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');
axis off equal; 
% caxis([min(10*floor(dep/10)),maxdep]);
caxis([minrep maxrep])
view(90,0);
material([0.5 0.5 0.5])
camlight left
if posi(3)>.5
	if (node>0),line(VER(node,1),VER(node,2),VER(node,3),'marker','o','markersize',8,'color','w','MarkerFaceColor','w','linestyle','none'); end
	if (onode>0),line(VER(onode,1),VER(onode,2),VER(onode,3),'marker','s','markersize',8,'color','w','MarkerFaceColor','w','linestyle','none'); end
	contourlines(VER,ITRI,rep,'delta',10,'labels',0)  
else
	contourlines(VER,ITRI,rep,'delta',20,'labels',0)  
end
axes('Position',[posi(1)+posi(3)*0.85 posi(2)+posi(4)*0.2 0.1*posi(3) 0.6*posi(4)]); axis off
caxis([minrep maxrep])
colorbar;

%%

