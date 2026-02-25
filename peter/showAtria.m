function showAtria(varargin)


node=0;onode=0;marker=0;
if length(varargin)<3
	error('Not enough parameters!!!');
else
	VER=varargin{1};
	ITRI=varargin{2};
	dep=varargin{3};
    endoVER=[];
	if (size(dep,1)==1),
        dep=dep';
    end	
	maxdep=max(dep);
	mindep=min(dep);
	maxi=maxdep;
	mini=mindep;
	lines=[];posi=[];
	pp=4;
	isolines=1;
	sublines=10;
    labels=0;    
	myview1=[90,0];
	myview=[-45,45];
    myview=[-63,15];
	sym=0;
	len=0;
	zoomfact=1;
	while pp<=length(varargin)
		if ischar(varargin{pp})
			key=lower(varargin{pp});
			switch key
				case 'nodes'
					node=varargin{pp+1};pp=pp+2;
				case 'marker'
					marker=varargin{pp+1};pp=pp+2;
				case 'isolines'
					isolines=varargin{pp+1};pp=pp+2;
				case 'labels'
					labels=varargin{pp+1};pp=pp+2;
				case 'onodes'
					onode=varargin{pp+1};pp=pp+2;
				case 'view'
					myview=varargin{pp+1};pp=pp+2;
				case 'view1'
					myview1=varargin{pp+1};pp=pp+2;					
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
                case 'endo'
                    endoVER=varargin{pp+1};pp=pp+2;                    
				case 'max'	
					maxi=varargin{pp+1};pp=pp+2;
					if length(maxi)==2
						mini=maxi(1);maxi=maxi(2);
                        givenMini=1;
					end				
				case 'range'										
					maxi=varargin{pp+1};pp=pp+2;
					if length(maxi)==2
						mini=maxi(1);maxi=maxi(2);
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
	if sym==0
% 		lines=loadmat('tims.mcm');
		lines=loadmat('hsvlim.mcm');
% 		lines=colormap('jet');		lines=lines(end:-1:1,:);
% 		A=colormap('hsv');			lines=A(1:45,:);colormap(lines);		
	else
		lines=loadmat('mymap_basic.mcm');	lines=lines(end:-1:1,:);
	end
end
if abs(mindep)>10  && ~exist('givenMini','var')
	mini=min(10*floor(mindep/10));
end
if sym && ~exist('givenMini','var')
	mini =-max(abs(dep));
	maxi = max(abs(dep));
end
if mini==maxi
    maxi=mini+1;
end
colormap(lines);
if isempty(posi)
	clf
	pos=get(gcf,'Position');
% 	set(gcf,'Position',[pos(1)  pos(2) 800   400]);
	posi=[0 0 1 1];
1end

if posi(3)>.5
	axes('Position',[posi(1)+posi(3)*0.42 posi(2)+posi(4)*0.2 0.1*posi(3) 0.6*posi(4)]); axis off
	caxis([mini maxi])
	colorbar;
end

%%
axes('Position',[posi(1)+0.01 posi(2) 0.44*posi(3) 1*posi(4)]);
patch('Faces',ITRI,'Vertices',VER,...
	  'FaceLighting','phong','BackFaceLighting','reverselit ','AmbientStrength',0.7,...
	  'FaceVertexCData',dep,'FaceColor','interp',...
	  'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');
axis off equal tight; 
caxis([mini maxi]);
view(myview1);
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
					 [VER(node(i),3) VER(node(i),3)-ci(3)*len],'marker','none','markersize',6,'color','k','MarkerFaceColor',[0.99 0.99 0.99],'linestyle','-','linewidth',2); 
			end
		else
			line(VER(node,1),VER(node,2),VER(node,3),'marker','o','markersize',6,'color','k','MarkerFaceColor',[0.99 0.99 0.99],'linestyle','none'); 
		end	
	end
	if (marker>0),line(VER(marker,1),VER(marker,2),VER(marker,3),'marker','o','markersize',15,'color',[0.99 0.99 0.99],'MarkerFaceColor',[0.99 0.99 0.99],'linestyle','none'); 	end
% 	if (onode>0),line(VER(onode,1),VER(onode,2),VER(onode,3),'marker','o','markersize',4,'color',[0.99 0.99 0.99],'MarkerFaceColor',[0.99 0.99 0.99],'linestyle','none'); end
	if (onode>0),line(VER(onode,1),VER(onode,2),VER(onode,3),'marker','o','markersize',10,'color',[0.99 0.99 0.99],'MarkerEdgeColor',[0.99 0.99 0.99],'linestyle','none'); end	
	if isolines,contourlines(VER,ITRI,dep,'delta',sublines,'labels',labels,'odd',1);  end
else
	if isolines,contourlines(VER,ITRI,dep,'delta',20,'labels',labels,'odd',1);  end
end

if ~isempty(endoVER)
    DrawLines(VER,ITRI,endoVER);
end
zoom(zoomfact);
camlight head
%%
axes('Position',[posi(1)+0.45*posi(3) posi(2) 0.6*posi(3) 1*posi(4)]);
patch('Faces',ITRI,'Vertices',VER,...
	  'FaceLighting','phong','BackFaceLighting','reverselit ','AmbientStrength',0.7,...
	  'FaceVertexCData',dep,'FaceColor','interp',...
	  'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');
axis off equal tight; 
caxis([mini maxi])
view(myview);

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
					 [VER(node(i),3) VER(node(i),3)-ci(3)*len],'marker','none','markersize',6,'color','k','MarkerFaceColor',[0.99 0.99 0.99],'linestyle','-','linewidth',2); 
			end
		else
			line(VER(node,1),VER(node,2),VER(node,3),'marker','o','markersize',6,'color','k','MarkerFaceColor',[0.99 0.99 0.99],'linestyle','none'); 
		end			
% 		line(VER(node,1),VER(node,2),VER(node,3),'marker','o','markersize',8,'color',[0.99 0.99 0.99],'MarkerFaceColor',[0.99 0.99 0.99],'linestyle','none'); 
	end
% 	if (onode>0),line(VER(onode,1),VER(onode,2),VER(onode,3),'marker','o','markersize',4,'color',[0.99 0.99 0.99],'MarkerFaceColor',[0.99 0.99 0.99],'linestyle','none'); end
if (marker>0),line(VER(marker,1),VER(marker,2),VER(marker,3),'marker','o','markersize',15,'color',[0.99 0.99 0.99],'MarkerFaceColor',[0.99 0.99 0.99],'linestyle','none'); 	end	
if (onode>0),line(VER(onode,1),VER(onode,2),VER(onode,3),'marker','s','markersize',8,'color',[0.99 0.99 0.99],'MarkerEdgeColor',[0.99 0.99 0.99],'linestyle','none'); end	
	if isolines,contourlines(VER,ITRI,dep,'delta',sublines,'labels',labels,'odd',1);  end
else
	if isolines, contourlines(VER,ITRI,dep,'delta',20,'labels',labels,'odd',1);  end
end
if ~isempty(endoVER)
    DrawLines(VER,ITRI,endoVER);
%     buur=graphdist(ITRI);
%     for i=1:length(buur)
%             buur(i,endoVER==endoVER(i))=0;
%     end
%     A=sum(buur)>0;
%     doITRI=zeros(length(ITRI),1);
%     for i=1:length(doITRI)
%         if sum(A(ITRI(i,:)))==3
%             doITRI(i)=1;
%         end
%     end
%    patch('Faces',ITRI(doITRI==1,:),'Vertices',VER,'Facecolor',[0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5])
end
zoom(zoomfact);
camlight head
%%
function DrawLines(VER,ITRI,endoVER)

    buur=graphdist(ITRI);
    adj=zeros(size(buur));
    for i=1:length(adj)
        endo=length(find(buur(i,endoVER==1)))>0;
        epi =length(find(buur(i,endoVER==0)))>0;
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