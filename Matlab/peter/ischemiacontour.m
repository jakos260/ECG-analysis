function ischemiacontour(varargin)


if length(varargin) < 2
	error('need the depolarization moments (1) and the adajcency matrix (2)');
else
	VER=varargin{1};
	ITRI=varargin{2};
	imi=[];
	if length(varargin) > 3
		pp=3;
		if ischar(varargin{pp})
			while pp<=length(varargin)
				if ischar(varargin{pp})
					key=lower(varargin{pp});
					switch key
						case 'imi'
							imi=varargin{pp+1};pp=pp+2;
						otherwise
							error('unknown parameter');
					end
				end				
			end
		elseif length(varargin{4})==1
			sym=varargin{4};
			if length(varargin) >4 && length(varargin{5})==1
				dolabels=varargin{5};
				if length(varargin) >5 && length(varargin{6})==1
					Nlevels=varargin{6};
				end
			end
			if length(varargin) >6
				disp('Paremeters 7 and further are ignored')
			end
		end		
	end
end

imiI=[];
for i=1:length(ITRI)
	if sum(ITRI(i,1)==imi | ITRI(i,2)==imi | ITRI(i,3)==imi)>2
		imiI=[imiI i];
	end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot contour lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx=14;
dz=20;
clf

axx=1:3;
myaxis=1;
maxax=max(max(abs(VER(:,axx(axx~=myaxis)))));
maxlev=range(VER(:,axx(axx==myaxis)));
for i=1:4
	axes('Position',[.25*(i-1) 0 0.25 0.5])
	level=maxlev(1)+25+(i-1)*dx;
	drawContour(VER,ITRI,VER(:,myaxis),level,'k',1);
	drawContour(VER,ITRI(imiI,:),VER(:,myaxis),level,'r',3);
	view(90,0);axis off;axis([-maxax maxax -maxax maxax -maxax maxax]);
end


myaxis=3;
maxax=max(max(abs(VER(:,axx(axx~=myaxis)))));
maxlev=range(VER(:,axx(axx==myaxis)));
for i=1:4
	axes('Position',[.25*(i-1) 0.5 0.25 0.5])
	level=maxlev(1)+15+(i-1)*dz;
	drawContour(VER,ITRI,VER(:,myaxis),level,'k',1);
	drawContour(VER,ITRI(imiI,:),VER(:,myaxis),level,'r',3);
	view(90,90);axis off;axis([-maxax maxax -maxax maxax -maxax maxax]);
end

%%

axes('Position',[.0 0.75 0.25 0.25])
patch('Faces',ITRI,'Vertices',VER,...
	  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
	  'FaceColor','g',...
	  'edgecolor','none','FaceAlpha',1);
patch('Faces',ITRI(imiI,:),'Vertices',VER,...
	  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
	  'FaceColor','r',...
	  'edgecolor','none','FaceAlpha',1);

  
axis off equal; 
material([0.5 0.5 0.5])
camlight left

myaxis=1;
maxax=max(max(abs(VER(:,axx(axx~=myaxis)))));
maxlev=range(VER(:,axx(axx==myaxis)));
for i=1:4
	level=maxlev(1)+25+(i-1)*dx;
	drawContour(VER,ITRI,VER(:,myaxis),level,'k',1);
end


myaxis=3;
maxax=max(max(abs(VER(:,axx(axx~=myaxis)))));
maxlev=range(VER(:,axx(axx==myaxis)));
for i=1:4
	level=maxlev(1)+15+(i-1)*dz;
	drawContour(VER,ITRI,VER(:,myaxis),level,'k',1);
end
view(90,0)
camlight left

axes('Position',[.25 0.75 0.25 0.25])
patch('Faces',ITRI,'Vertices',VER,...
	  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
	  'FaceColor','g',...
	  'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');
material([0.5 0.5 0.5])
patch('Faces',ITRI(imiI,:),'Vertices',VER,...
	  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
	  'FaceColor','r',...
	  'edgecolor','none','FaceAlpha',1);
axis off equal; 



myaxis=1;
maxax=max(max(abs(VER(:,axx(axx~=myaxis)))));
maxlev=range(VER(:,axx(axx==myaxis)));
for i=1:4
	level=maxlev(1)+25+(i-1)*dx;
	drawContour(VER,ITRI,VER(:,myaxis),level,'k',1);
end


myaxis=3;
maxax=max(max(abs(VER(:,axx(axx~=myaxis)))));
maxlev=range(VER(:,axx(axx==myaxis)));
for i=1:4
	level=maxlev(1)+15+(i-1)*dz;
	drawContour(VER,ITRI,VER(:,myaxis),level,'k',1);
end
view(180,0)
camlight left

%%------------------------------------------------------
function drawContour(VER,ITRI,fun,levels,contourcolor,width)

triangle_val = fun(ITRI);
triangle_min = min(triangle_val, [], 2);
triangle_max = max(triangle_val, [], 2);

for cnt_indx=1:length(levels)
	cnt = levels(cnt_indx);
	use = cnt>=triangle_min & cnt<=triangle_max;
	intersect1 = [];
	intersect2 = [];
	
	for tri_indx=find(use)'
		pos  = VER(ITRI(tri_indx,:), :);
		v(1) = triangle_val(tri_indx,1);
		v(2) = triangle_val(tri_indx,2);
		v(3) = triangle_val(tri_indx,3);
		la(1) = (cnt-v(1)) / (v(2)-v(1)+eps); % abcissa between vertex 1 and 2
		la(2) = (cnt-v(2)) / (v(3)-v(2)+eps); % abcissa between vertex 2 and 3
		la(3) = (cnt-v(3)) / (v(1)-v(3)+eps); % abcissa between vertex 1 and 2
		
		b=cross(pos(2,:) - pos(1,:),pos(3,:) - pos(1,:));
		c=mean(b); c=eps*c./norm(c);
		abc(1,:) = pos(1,:) + la(1) * (pos(2,:) - pos(1,:)) + c;
		abc(2,:) = pos(2,:) + la(2) * (pos(3,:) - pos(2,:)) + c;
		abc(3,:) = pos(3,:) + la(3) * (pos(1,:) - pos(3,:)) + c;
		sel     = find(la>=0 & la<=1);
		if size(sel,2)>=2,
			intersect1 = [intersect1; abc(sel(1),:)];
			intersect2 = [intersect2; abc(sel(2),:)];
		end
	end

	% remember the details for external reference
	contour(cnt_indx).level = cnt;
	contour(cnt_indx).n     = size(intersect1,1);
	contour(cnt_indx).intersect1 = intersect1;
	contour(cnt_indx).intersect2 = intersect2;

end

% collect all different contourlevels and plot them
intersect1 = [];
intersect2 = [];
cntlevel   = [];

for cnt_indx=1:length(levels)
  intersect1 = [intersect1; contour(cnt_indx).intersect1];
  intersect2 = [intersect2; contour(cnt_indx).intersect2];
  cntlevel   = [cntlevel; ones(contour(cnt_indx).n,1) * levels(cnt_indx)];

end

if  isempty(intersect1)|| isempty(intersect2)  
	return;
end
X = [intersect1(:,1) intersect2(:,1)]';
Y = [intersect1(:,2) intersect2(:,2)]';
C = [cntlevel(:)     cntlevel(:)]';

if size(VER,2)>2
  Z = [intersect1(:,3) intersect2(:,3)]';
else
  Z = zeros(2, length(cntlevel));
end
white=[ 0.99 0.99 0.99];
hc = [];
for i=1:length(cntlevel)
	h1 = line(X(:,i), Y(:,i), Z(:,i),'color',contourcolor,'userdata',cntlevel(i),'linestyle','-','linewidth',width);
    hc = [hc; h1];
end

%%
function drawIschemiaContour(VER,ITRI,fun,levels,contourcolor,width)

trans=1;
triangle_val = fun(ITRI);
triangle_min = min(triangle_val, [], 2);
triangle_max = max(triangle_val, [], 2);

for cnt_indx=1:length(levels)
	cnt = levels(cnt_indx);
	use = cnt>=triangle_min & cnt<=triangle_max;
	intersect1 = [];
	intersect2 = [];
	
	for tri_indx=find(use)'
		pos  = VER(ITRI(tri_indx,:), :);
		v(1) = triangle_val(tri_indx,1);
		v(2) = triangle_val(tri_indx,2);
		v(3) = triangle_val(tri_indx,3);
		la(1) = (cnt-v(1)) / (v(2)-v(1)+eps); % abcissa between vertex 1 and 2
		la(2) = (cnt-v(2)) / (v(3)-v(2)+eps); % abcissa between vertex 2 and 3
		la(3) = (cnt-v(3)) / (v(1)-v(3)+eps); % abcissa between vertex 1 and 2
		
		b=cross(pos(2,:) - pos(1,:),pos(3,:) - pos(1,:));
		c=mean(b); c=eps*c./norm(c);
		abc(1,:) = pos(1,:) + la(1) * (pos(2,:) - pos(1,:)) + c;
		abc(2,:) = pos(2,:) + la(2) * (pos(3,:) - pos(2,:)) + c;
		abc(3,:) = pos(3,:) + la(3) * (pos(1,:) - pos(3,:)) + c;
		sel     = find(la>=0 & la<=1);
		if size(sel,2)>=2,
			intersect1 = [intersect1; abc(sel(1),:)];
			intersect2 = [intersect2; abc(sel(2),:)];
		end
	end

	% remember the details for external reference
	contour(cnt_indx).level = cnt;
	contour(cnt_indx).n     = size(intersect1,1);
	contour(cnt_indx).intersect1 = intersect1;
	contour(cnt_indx).intersect2 = intersect2;

end

% collect all different contourlevels and plot them
intersect1 = [];
intersect2 = [];
cntlevel   = [];

for cnt_indx=1:length(levels)
  intersect1 = [intersect1; contour(cnt_indx).intersect1];
  intersect2 = [intersect2; contour(cnt_indx).intersect2];
  cntlevel   = [cntlevel; ones(contour(cnt_indx).n,1) * levels(cnt_indx)];

end

if  isempty(intersect1)|| isempty(intersect2)  
	return;
end
X = [intersect1(:,1) intersect2(:,1)]';
Y = [intersect1(:,2) intersect2(:,2)]';
C = [cntlevel(:)     cntlevel(:)]';

if size(VER,2)>2
  Z = [intersect1(:,3) intersect2(:,3)]';
else
  Z = zeros(2, length(cntlevel));
end
white=[ 0.99 0.99 0.99];
hc = [];
for i=1:length(cntlevel)
	h1 = line(X(:,i), Y(:,i), Z(:,i),'color',contourcolor,'userdata',cntlevel(i),'linestyle','-','linewidth',width);

    hc = [hc; h1];
end


