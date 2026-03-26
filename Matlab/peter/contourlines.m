function contourlines(varargin)% VER,ITRI,VALS,sym)


if length(varargin) < 3
	error('need the depolarization moments (1) and the adajcency matrix (2)');
else
	VER=varargin{1};
	ITRI=varargin{2};
	VALS=varargin{3};
	dolabels=0;
	sym=0;
	Nlevels=-1;
	delta=Inf;
	extremes=0;
	scale='lin';
	doodd=0;
	tcolor=[0.99 0.99 0.99];
	if length(varargin) > 3
		pp=4;
		if ischar(varargin{pp})
			while pp<=length(varargin)
				if ischar(varargin{pp})
					key=lower(varargin{pp});
					switch key
						case 'delta'
							delta=varargin{pp+1};pp=pp+2;
						case 'sym'
							sym=varargin{pp+1};pp=pp+2;
						case 'labels'
							dolabels=varargin{pp+1};pp=pp+2;
						case 'extremes'
							extremes=varargin{pp+1};pp=pp+2;
						case 'color'
							tcolor=varargin{pp+1};pp=pp+2;
						case 'scale'
							scale=varargin{pp+1};pp=pp+2;
% 							if ~foundlabels, dolabels=0;end;
						case 'odd'
							doodd=varargin{pp+1};pp=pp+2;
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


if min(VALS) ==max(VALS)
    return;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot contour lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contourcolor='k';
fun=VALS;
triangle_val = fun(ITRI);
triangle_min = min(triangle_val, [], 2);
triangle_max = max(triangle_val, [], 2);
absmax = max(abs([min(VALS) max(VALS)]));mabsmax=ceil(log10(absmax))-1;
maxi=max(VALS);if maxi,maxim=real(ceil(log10(maxi))-1);		 else maxim=-10000; end
mini=min(VALS);if mini, minim=real(ceil(log10(abs(mini)))-1);else minim=-10000; end
mantisse=max(minim,maxim);
tmaxi=floor(maxi/10^mantisse);
tmini=floor(mini/10^mantisse);
range=tmaxi-tmini;
if sym
    tmaxi=floor(absmax/10^mabsmax);
    tmini=-floor(absmax/10^mabsmax);
end
if delta~=Inf
	zebra=delta;
	levels= ceil(mini/zebra)*zebra:zebra:floor(maxi/zebra)*zebra;
elseif Nlevels>0
    zebra=(maxi-mini)/Nlevels;
	levels = mini+zebra:zebra:maxi-zebra;
elseif strcmp(scale,'log')
	levels=[];j=1;
	for i=min(minim,maxim-1):1:maxim
		levels(j)=10^i;
		j=j+1;
	end
	levels=[-levels(end:-1:1) 0 levels];	
else
	
	zebra=0.1*ceil(range/10)*10^mantisse;
	zebraOrg=zebra;
	if tmini*10^mantisse <= 0 && tmaxi*10^mantisse >=0
		aaa=0:-zebra:tmini*10^mantisse;
		levels=[aaa(end:-1:1) 0:zebra:tmaxi*10^mantisse];
	else
		levels = tmini*10^mantisse:zebra:tmaxi*10^mantisse;
	end
	bbb=[1 2 2.5 5 10 20 25 50 100 200 250 500 1000];
% 	bbb=[1 10 100 1000];	
	j=1;
	while length(levels) >10 && j<=length(bbb)
		zebra=zebraOrg*bbb(j);
		j=j+1;
		if tmini*10^mantisse <= 0 && tmaxi*10^mantisse >=0
			aaa=0:-zebra:tmini*10^mantisse;
			levels=[aaa(end:-1:1) zebra:zebra:tmaxi*10^mantisse];
		else
			levels = tmini*10^mantisse:zebra:tmaxi*10^mantisse;
		end
    end
    levels(levels+0.5*zebra>maxi)=[];
    levels(levels-0.5*zebra<mini)=[];
end
if max(VALS) == maxi
    levels(levels==maxi)=[];
end

if doodd %length(VALS(VALS==0))> length(VALS)/3
	levels(levels==0)=[];
    
end



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
	if cntlevel(i)>0
	    h1 = line(X(:,i), Y(:,i), Z(:,i),'color',contourcolor,'userdata',cntlevel(i),'linestyle','-');
	elseif cntlevel(i)<0
		h1 = line(X(:,i), Y(:,i), Z(:,i),'color',contourcolor,'userdata',cntlevel(i),'linestyle','--');
	else
		h1 = line(X(:,i), Y(:,i), Z(:,i),'color',contourcolor,'userdata',cntlevel(i),'linestyle',':','linewidth',1);
	end
    hc = [hc; h1];
end



if dolabels
	i=1;	
	deltalab=20;
% 	ib=cntlevel(i);	
	iib=i;
	ib=-999999;
	deltal=deltalab;
	while i<length(cntlevel) %&& cntlevel(i)~=cntlevel(end)
		deltal=deltal+1;
		a=find(cntlevel==cntlevel(i));
		if length(a) >10 && ( deltal>=deltalab)%&& cntlevel(i)~=cntlevel(end)%|| i==length(cntlevel)
			deltal=0;
			j=round((i+iib)/2);
			
			x=mean(X(1,j));    y=mean(Y(1,j));	z=mean(Z(1,j));
			B=norm3d([VER(:,1)-x VER(:,2)-y VER(:,3)-z]);
			A=find(B==min(B));A=A(1);
			[ti,l]=find(ITRI==A);
			b=cross(VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),2),:),...
				VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),3),:));
			c=mean(b); c=c./norm(c);
			c=get(gca,'CameraTarget')-get(gca,'CameraPosition');c=-70*c/norm(c);
			x=x+c(1);y=y+c(2);z=z+c(3);
			
			if abs(cntlevel(iib)) >1e-15,
				txt=num2str(cntlevel(iib));
			else
				txt='0';
			end
			text(x,y,z,txt,'HorizontalAlignment','center','EraseMode','normal','color',tcolor,'fontname','verdana','fontsize',6);
			a=find(cntlevel==cntlevel(i));
			i=a(end)+1;
			if i<length(cntlevel)
				a=find(cntlevel==cntlevel(i));
				i=a(end)+1;
			end		
			if i<length(cntlevel)
				a=find(cntlevel==cntlevel(i));
				ib=cntlevel(i);
				iib=i;

				if (length(a) >5  )%&& ib~=cntlevel(end))|| (length(a) >25  && ib==cntlevel(end))
					deltalab=length(a)/3;
					deltal=deltalab;
				else
					deltalab=100000;
				end
			end
		end
		i=i+1;
	end
end

% extreme values
if extremes
	maxVer=find(VALS==max(VALS));maxVer=maxVer(1);maxVal=max(VALS);
	minVer=find(VALS==min(VALS));minVer=minVer(1);minVal=min(VALS);
	if 1%length(maxVer)==1
		[ti,l]=find(ITRI==maxVer);
		b=cross(VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),2),:),...
		VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),3),:));
		c=mean(b); c=-5*c./norm(c);
		c=get(gca,'CameraTarget')-get(gca,'CameraPosition');c=-71*c/norm(c);		
		x=VER(maxVer,1)+c(1);y=VER(maxVer,2)+c(2);z=VER(maxVer,3)+c(3);
			
		if abs(maxVal) >1e-15,
			if abs(maxVal) < 1 || (abs(maxVal) <100 && abs(maxVal)>10)
				txt=num2str(maxVal,2);
			else
				txt=num2str(maxVal,3);
			end
			if maxVal >0, txt=['+' txt];end
		else
			txt='0';
		end
		text(x,y,z,txt,'HorizontalAlignment','left','EraseMode','normal','Fontweight','bold','Color',tcolor,'fontname','verdana','fontsize',7);
	end
	if 1%length(minVer)==1
		[ti,l]=find(ITRI==minVer);
		b=cross(VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),2),:),...
		VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),3),:));
		c=mean(b); c=-7*c./norm(c);
		c=get(gca,'CameraTarget')-get(gca,'CameraPosition');c=-71*c/norm(c);
		
		x=VER(minVer,1)+c(1);y=VER(minVer,2)+c(2);z=VER(minVer,3)+c(3);
			
		if abs(minVal) >1e-15,
			if abs(minVal) <1 || (abs(minVal) <100 && abs(minVal)>10)
				txt=num2str(minVal,2);
			else
				txt=num2str(minVal,3);
			end
		else
			txt='0';
		end
		text(x,y,z,txt,'HorizontalAlignment','left','EraseMode','normal','Fontweight','bold', 'Color',tcolor,'fontname','verdana','fontsize',7);
	end
% text(VER(maxVer,1),VER(maxVer,2),VER(maxVer,3),['+' num2str(maxVal,3)],'HorizontalAlignment','center')
% text(VER(minVer,1),VER(minVer,2),VER(minVer,3),['-' num2str(minVal,3)],'HorizontalAlignment','center');
end