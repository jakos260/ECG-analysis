function distline(varargin)% VER,ITRI,VALS,sym)


if length(varargin) < 3
	error('need the depolarization moments (1) and the adajcency matrix (2)');
else
	VER=varargin{1};
	ITRI=varargin{2};
	REF=varargin{3};
	levels=varargin{4};
	adj=varargin{5};
	VALS=graphdistone(adj,REF);
end


dolabels=1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot contour lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contourcolor='y';
fun=VALS;
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

hc = [];
for i=1:length(cntlevel)
    h1 = patch('XData', X(:,i), 'Ydata', Y(:,i), ...
               'ZData', Z(:,i), 'CData', C(:,i), ...
               'facecolor','none','edgecolor',contourcolor,...
               'userdata',cntlevel(i));
    hc = [hc; h1];
end
if dolabels
	i=1;	ib=cntlevel(i);	iib=1;
	while i<=length(cntlevel)
		if cntlevel(i)~= ib || i==length(cntlevel)
			j=round((i+iib)/2);
			x=mean(X(:,j));    y=mean(Y(:,j));	z=mean(Z(:,j));
			B=norm3d([VER(:,1)-x VER(:,2)-y VER(:,3)-z]);
			A=find(B==min(B));
			[ti,l]=find(ITRI==A);
			b=cross(VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),2),:),...
				VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),3),:));
			c=mean(b); c=-5*c./norm(c);
			x=x+c(1);y=y+c(2);z=z+c(3);
			
			if abs(cntlevel(iib)) >1e-15,
				txt=num2str(cntlevel(iib));
			else
				txt='0';
			end
% 			text(x,y,z,txt,'HorizontalAlignment','center','EraseMode','normal','color',contourcolor);
			a=find(cntlevel==cntlevel(i));
			i=a(end)+1;
			if i<length(cntlevel)
				ib=cntlevel(i);
				iib=i;
			end
		end
		i=i+1;
	end
end

[ti,l]=find(ITRI==REF);
b=cross(VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),2),:),...
VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),3),:));
c=mean(b); c=-5*c./norm(c);
c=VER(REF,:)+c;
line(c(1),c(2),c(3),'Marker','o','color','w','markerfacecolor','w')