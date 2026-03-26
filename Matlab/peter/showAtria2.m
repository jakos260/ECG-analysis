function showAtria2(varargin)


node=0;onode=0;
if length(varargin)<3
	error('Not enough parameters!!!');
else
	VER=varargin{1};
	ITRI=varargin{2};
	dep1=varargin{3};
	dep2=varargin{4};
    endoVER=[];
	if (size(dep1,1)==1),dep1=dep1';end	
	if (size(dep2,1)==1),dep2=dep2';end		
	deps=[dep1; dep2];
	maxdep=max(deps);
	mindep=min(deps);
	maxi=maxdep;
	mini=mindep;
	lines=[];posi=[];myview=[];
	pp=5;
	isolines=1;
	sublines=10;
    labels=0;    
	myview=[90,0];
	sym=0;
	do2=0;
	len=0;
	zoomfact=1;
	while pp<=length(varargin)
		if ischar(varargin{pp})
			key=lower(varargin{pp});
			switch key
				case 'nodes'
					node=varargin{pp+1};pp=pp+2;
				case 'isolines'
					isolines=varargin{pp+1};pp=pp+2;
				case 'labels'
					labels=varargin{pp+1};pp=pp+2;
				case 'onodes'
					onode=varargin{pp+1};pp=pp+2;
				case 'view'
					myview=varargin{pp+1};pp=pp+2;
                case 'endo'
                    endoVER=varargin{pp+1};pp=pp+2;                					
				case 'sublines'
					sublines=varargin{pp+1};pp=pp+2;
				case 'sym'
					sym=varargin{pp+1};pp=pp+2;
				case 'colormap'
					lines=varargin{pp+1};pp=pp+2;					
				case 'tail'
					len=varargin{pp+1};pp=pp+2;					
				case 'zoom'					
					zoomfact=varargin{pp+1};pp=pp+2;
				case 'do2'					
					do2=varargin{pp+1};pp=pp+2;					
				case 'max'					
					maxi=varargin{pp+1};pp=pp+2;
					if length(maxi)==2
						mini=maxi(1);maxi=maxi(2);
                        givenMini=1;
					elseif length(maxi)==4
						mini=maxi([1,3]);
						maxi=maxi([2,4]);
                        givenMini=1;						
					end
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
set(gcf,'PaperPositionMode','auto')

if isempty(lines)
% 	lines=loadmat('tims.mcm');
% 	lines=colormap('hsv');		lines=lines(end:-1:1,:);
% 	A=colormap('hsv');	lines=A(1:45,:);colormap(lines);
	lines=loadmat('hsvlim.mcm');%lines=lines(end:-1:1,:);
end
if abs(mindep)>10  && ~exist('givenMini')
	mini=min(10*floor(mindep/10));
end
if sym && ~exist('givenMini')
	mini=-max(abs(deps));
	maxi=max(abs(deps));
end

colormap(lines);
if isempty(posi)
	clf
% 	pos=get(gcf,'Position');
% 	set(gcf,'Position',[pos(1)  pos(2) 800   400]);
	posi=[0 0 1 1];
end

if posi(3)>.5 && ~do2
	axes('Position',[posi(1)+posi(3)*0.42 posi(2)+posi(4)*0.2 0.1*posi(3) 0.6*posi(4)]); axis off
	caxis([mini(1) maxi(1)])
	colorbar;
end

%%
axes('Position',[posi(1)+.01 posi(2) 0.45*posi(3) 1*posi(4)]);
patch('Faces',ITRI,'Vertices',VER,...
	  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
	  'FaceVertexCData',dep1,'FaceColor','interp',...
	  'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');
if ~isempty(endoVER)
    DrawLines(VER,ITRI,endoVER);
end

  axis off equal tight; 
if do2
	if exist('givenMini') && givenMini
		caxis([mini(1) maxi(1)])
	else
		caxis([min(dep1) max(dep1)])
	end
else
	caxis([mini(1) maxi(1)])	
end
view(myview);
zoom(zoomfact);
material([0.5 0.5 0.5])
% camlight left
if posi(3)>.5
	if (node>0),
		if len>0
			for i=1:length(node)
				[ti,l]=find(ITRI==node(i));
				b=cross(VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),2),:),VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),3),:));
				c=mean(b); ci=c./norm(c); % mm
				line([VER(node(i),1) VER(node(i),1)-ci(1)*len],...
					 [VER(node(i),2) VER(node(i),2)-ci(2)*len],...
					 [VER(node(i),3) VER(node(i),3)-ci(3)*len],'marker','none','markersize',8,'color','k','MarkerFaceColor',[0.99 0.99 0.99],'linestyle','-','linewidth',2); 
			end
		else
			line(VER(node,1),VER(node,2),VER(node,3),'marker','o','markersize',8,'color','k','MarkerFaceColor',[0.99 0.99 0.99],'linestyle','none'); 
		end	
	end
	if (onode>0),line(VER(onode,1),VER(onode,2),VER(onode,3),'marker','o','markersize',6,'color',[0.99 0.99 0.99],'MarkerFaceColor',[0.01 0.01 0.01],'linestyle','none'); end
	if isolines,contourlines(VER,ITRI,dep1,'delta',sublines,'labels',labels,'odd',1);  end
else
	if isolines,contourlines(VER,ITRI,dep1,'delta',20,'labels',labels,'odd',1);  end
end
if do2,
	colorbar('location','southoutside');
end

%%
axes('Position',[posi(1)+0.55*posi(3)-0.01 posi(2) 0.45*posi(3) 1*posi(4)-0.01]);
patch('Faces',ITRI,'Vertices',VER,...
	  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
	  'FaceVertexCData',dep2,'FaceColor','interp',...
	  'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');
if ~isempty(endoVER)
    DrawLines(VER,ITRI,endoVER);
end
axis off equal tight; 
if do2
	if exist('givenMini') && givenMini
		caxis([mini(end) maxi(end)])
	else
		caxis([min(dep1) max(dep1)])
	end
else
	caxis([mini(end) maxi(end)])	
end
view(myview);
zoom(zoomfact);
material([0.5 0.5 0.5])
% camlight left
if posi(3)>.5
	if (node>0),
	if len>0
			for i=1:length(node)
				[ti,l]=find(ITRI==node(i));
				b=cross(VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),2),:),VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),3),:));
				c=mean(b); ci=c./norm(c); % mm
				line([VER(node(i),1) VER(node(i),1)-ci(1)*len],...
					 [VER(node(i),2) VER(node(i),2)-ci(2)*len],...
					 [VER(node(i),3) VER(node(i),3)-ci(3)*len],'marker','none','markersize',8,'color','k','MarkerFaceColor',[0.99 0.99 0.99],'linestyle','-','linewidth',2); 
			end
		else
			line(VER(node,1),VER(node,2),VER(node,3),'marker','o','markersize',8,'color','k','MarkerFaceColor',[0.99 0.99 0.99],'linestyle','none'); 
		end			
% 		line(VER(node,1),VER(node,2),VER(node,3),'marker','o','markersize',8,'color',[0.99 0.99 0.99],'MarkerFaceColor',[0.99 0.99 0.99],'linestyle','none'); 
	end
	if (onode>0),line(VER(onode,1),VER(onode,2),VER(onode,3),'marker','o','markersize',6,'color',[0.99 0.99 0.99],'MarkerFaceColor',[0.01 0.01 0.01],'linestyle','none'); end
	if isolines,contourlines(VER,ITRI,dep2,'delta',sublines,'labels',labels,'odd',1);  end
else
	if isolines, contourlines(VER,ITRI,dep2,'delta',20,'labels',labels,'odd',1);  end
end
if do2,	colorbar('location','southoutside');end

%%
function DrawLines(VER,ITRI,endoVER)

    buur=graphdist(ITRI);
    adj=zeros(size(buur));
    for i=1:length(adj)
        endo=length(find(buur(i,endoVER==1)))>0;
        epi =length(find(buur(i,endoVER==0)))>0;
        if i==438
            stop=1;
        end
        if endo && epi && endoVER(i)==1
            a=find(buur(i,:)>0);
            a(endoVER(a)==0)=[];
            adj(i,a)=1;
            adj(a,i)=1;                   
        end
    end
    for i=1:length(adj)
        endo=length(find(buur(i,endoVER==1)))>0;
        epi =length(find(buur(i,endoVER==0)))>0;
        if endo ~=epi 
            adj(i,:)=0;
            adj(:,i)=0;                               
        end
    end
    A=splitgraph(adj);
    B=zeros(size(A));
    k=0;
    for i=1:max(A)
        if length(find(A==i))>1
            k=k+1;
            B(A==i)=k;
            A(A==i)=0;
        end
    end
    for i=1:k
        L=find(B==i);
        l0=L(1);l=l0;
        while length(L)>1
            r=getroute(adj,l,L(2));
            line(VER(r,1),VER(r,2),VER(r,3),'color',[0.99 0.99 0.99],'linewidth',1);
            l=r(end);
            for j=1:length(r)-1,L(L==r(j))=[];end
        end
        r=getroute(adj,l0,l);
        line(VER(r,1),VER(r,2),VER(r,3),'color',[0.99 0.99 0.99],'linewidth',1);
	end
	