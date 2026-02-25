 function showpatch(varargin)


node=0;onode=0;marker=0;
if length(varargin)<3
	error('Not enough parameters!!!');
else
	VER=varargin{1};
	ITRI=varargin{2};
    pp=3;
    alpha = 1;
    dep=[];
    if ~ischar(varargin{3})
        dep=varargin{3};
        pp=4;
    end
	endoVER=[];
	if (size(dep,1)==1),dep=dep';end	
	maxdep=max(dep);
	mindep=min(dep);
	lines=[];posi=[];
    myview=[];
	isolines=1;
	sublines=10;
	edges='none';
	dobar=0;
	doclf=0;
	myview=[90,0];
	color=[];
    elec =[];
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
				case 'onodes'
					onode=varargin{pp+1};pp=pp+2;
				case 'view'
					myview=varargin{pp+1};pp=pp+2;
				case 'sublines'
					sublines=varargin{pp+1};pp=pp+2;
				case 'edge'
					edges=varargin{pp+1};pp=pp+2;
				case 'dobar'
					dobar=varargin{pp+1};pp=pp+2;                   
				case 'alpha'
					alpha=varargin{pp+1};pp=pp+2;
				case 'range'
					drang=varargin{pp+1};pp=pp+2;
					mindep=drang(1);
					maxdep=drang(2);
				case 'clf'
					doclf=varargin{pp+1};pp=pp+2;
                case 'endo'
                    endoVER=varargin{pp+1};pp=pp+2;    	
				case 'position'	
					posi=varargin{pp+1};pp=pp+2; 
				case 'color'	
					color=varargin{pp+1};pp=pp+2; 
				case 'elec'	
					elec=varargin{pp+1};pp=pp+2; 
                    
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
% set(gcf,'PaperPositionMode','auto')
if isempty(lines)
	lines=loadmat('tims.mcm');
% 	A=colormap('hsv');	lines=A(1:45,:);colormap(lines);
% 	lines=loadmat('hsvlim.mcm');
% 	lines=colormap('jet');			lines=lines(end:-1:1,:);
end
if ~isempty(color)
	dobar=0;
	isolines=0;
end
if abs(mindep)>10
	mindep=min(10*floor(mindep/10));
end
colormap(lines);
if doclf
    clf; 
end
if isempty(posi)
	posi=[0 0 1 1];
end
if dobar
	axes('Position',[posi(1)+0.01 posi(2)+0.01 0.84*posi(3) 0.98*posi(4)]);
else
	axes('Position',[posi(1)+0.01 posi(2)+0.01 0.95*posi(3) 0.95*posi(4)]);
end

hold on;
if posi(3)>.5 && dobar
% 	axes('Position',[posi(1)+posi(3)*0.7 posi(2)+posi(4)*0.1 0.25*posi(3) 0.8*posi(4)]); axis off
	caxis([mindep maxdep])
	colorbar;
end

%%

if isempty(color)
	patch('Faces',ITRI,'Vertices',VER,...
		  'FaceLighting','phong','BackFaceLighting','reverselit','AmbientStrength',0.7,...
		  'FaceVertexCData',dep,'FaceColor','interp',...
		  'edgecolor',edges,'FaceAlpha',alpha,'buttondownFcn','selectnode');
	  caxis([mindep,maxdep]);
else
	patch('Faces',ITRI,'Vertices',VER,...
		  'FaceLighting','phong','BackFaceLighting','reverselit','AmbientStrength',0.7,...
		  'FaceColor',color,...
		  'edgecolor',edges,'FaceAlpha',alpha,'buttondownFcn','selectnode','parent',gca);
end
a= range(VER);a=a(:)' *1.3;
if a(1) ==a(2)
    a(2) = a(1)+1;
end
axis(a)

axis off vis3d;

view(myview);
material([0.5 0.5 0.5])
camlight left
if ~isempty(elec)
line(elec(:,1),elec(:,2),elec(:,3),'marker','o','markersize',10,'color',[0 0 0],'MarkerFaceColor',[0 0 0],'linestyle','none','parent',gca); 
end

if posi(3)>.5
	if (node>0),line(VER(node,1),VER(node,2),VER(node,3),'marker','o','markersize',8,'color','w','MarkerFaceColor','w','linestyle','none'); end
	if (onode>0),line(VER(onode,1),VER(onode,2),VER(onode,3),'marker','s','markersize',4,'color','k','MarkerFaceColor','k','linestyle','none'); end
	if (marker>0),line(VER(marker,1),VER(marker,2),VER(marker,3),'marker','o','markersize',16,'color',[0.99 0.99 0.99],'MarkerFaceColor',[0.99 0.99 0.99],'linestyle','none'); 	end		
	if isolines,contourlines(VER,ITRI,dep,'delta',sublines,'labels',0,'odd',1);  end
else
	if isolines,contourlines(VER,ITRI,dep,'delta',20,'labels',0,'odd',1);  end
end
if ~isempty(endoVER)
    DrawLines(VER,ITRI,endoVER);
end

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

%%
