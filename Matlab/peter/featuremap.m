function featuremap(varargin)
% clf
dmini=cast(intmin('int16'),'double');
if length(varargin) < 1
	error('This routine needs at least two parameters');
else
	sig=varargin{1};
	level=-1;
	maxi=max(max(abs(sig(sig~=dmini))));
	pp=2;
	markers=[];
	episodes=[];
	labels=[];
	dodiff=0;
	dozmap=0;
	deltaError=0;
	doline=1;
	while pp<=nargin
		if ischar(varargin{pp})
			key=lower(varargin{pp});
			switch key
				case 'level'
					level=varargin{pp+1};pp=pp+2;
					if ~isempty(level)
						maxi=2*level;
					else
						level=-1;
					end
				case 'max'
					maxi=varargin{pp+1};pp=pp+2;
				case 'labels'
					labels=varargin{pp+1};pp=pp+2;
				case 'markers'
					markers=varargin{pp+1};pp=pp+2;
				case 'episodes'
					episodes=varargin{pp+1};pp=pp+2;
				case 'dodiff'
					dodiff=varargin{pp+1};pp=pp+2;
				case 'dozmap'
					dozmap=varargin{pp+1};pp=pp+2;
				case 'doline'
					doline=varargin{pp+1};pp=pp+2;
				case 'xxx'
					deltaError=varargin{pp+1};pp=pp+2;
				otherwise
					error('unknown parameter');
			end
		end
	end
end
if dozmap
	org=sig;
	for i=1:size(sig,2)
		tmp=sig(:,i);tmp(isnan(tmp))=[];
		sig(:,i)=(sig(:,i)-mean(tmp))./std(tmp);
	end
	sig(isnan(org))=NaN;
	maxi=3;
end
if dodiff
	tmp=sig;
	a1=ones(size(sig,1),1)*sig(1,:);
	a1(a1==dmini)=0;	
% 	tmp(tmp==dmini)=0;		
	a1=(tmp-a1);%./a1;
	a1(sig==dmini)=dmini;
	a1(a1==-Inf)=NaN;
	a1(a1==Inf)=NaN;	
	sig=a1;
end

if deltaError~=0
	A=sig';
	iNAN=isnan(A);
	A(iNAN)=2*max(max(A));
	for i=1:rank(A)
		B=tsvd(A,i);
		if 100*max(max(abs(B-A)))<deltaError
			break;
		end
	end	
	B(iNAN)=NaN;
	sig=B';
end


if ~isempty(markers)
	x=markers;	
	if size(x,1)~=1,x=x';end
elseif ~isempty(episodes)
	x=1:length(episodes);
	if size(x,1)~=1,x=x';end	
else
	x=(0:size(sig,1)-1)*10;
end
y=1:size(sig,2);
yt=y;
if isempty(labels)
	labels=num2str(y');
end
if size(y,1)==1 && size(y,2)==1
	y=[1 2 3];
	sig=[sig,sig,sig];
	yt=2;
	is12lead=1;
else
	is12lead=0;
end



yg=y;
xg=x;
gsig=sig;
csig=sig;
reverseax=0;
if x(1) > x(end)
	reverseax=1;
end
if ~isempty(episodes)
	if x(1) > x(end)
		reverseax=1;
		x=sort([x x-length(x)/10000],'descend');
	else
		x=sort([x x+length(x)/10000]);
	end
	csig=zeros(length(x),length(y));
	k=1;
	csig(1,:)=sig(k,:);
	for i=2:2:size(csig,1)-1
		csig(i,:)=sig(k,:);
		csig(i+1,:)=sig(k,:);
		k=k+1;
	end
	csig(end,:)=sig(end,:);
	tmp=csig(:,1);
	ltmp=labels(1,:);
	sp=[];for i=1:length(ltmp),sp=[sp ' '];end
	for i=2:size(csig,2)
		tmp=[tmp, (csig(:,i-1)+csig(:,i))./2];
		tmp=[tmp, csig(:,i)];		
	end
	for i=2:length(yt)
		ltmp=[ltmp; sp];
		ltmp=[ltmp; labels(i,:)];		
	end
	labels=ltmp;
	csig=tmp;
	y=1:0.5:size(sig,2);
end

XG=xg'*ones(size(yg));
YG=(yg'*ones(size(xg)))';
X=x'*ones(size(y));
Y=(y'*ones(size(x)))';

if dozmap
	cols=loadmat('sigmamap.mcm');
else
	cols=loadmat('mymap_basic.mcm');
end
colormap(cols);
csig(csig==dmini)=NaN;

if (~isempty(episodes) )
	[ITRI,VER,dat]=surf2patch(X,Y,zeros(size(Y)),csig);
% 	dat(isnan(dat))=0;
	patch('Faces',ITRI,'Vertices',VER,'FaceVertexCData',dat,'FaceColor','interp','edgecolor','none','Parent',gca);
	
	hold on
	surface(XG,YG,zeros(size(YG)),gsig,'EdgeColor','k','FaceColor','none','FaceAlpha',0.1,'Parent',gca);
	hold off
	axis([min(x) max(x)  min(y) max(y) ])
else
	surface(X,Y,csig,'EdgeColor','none','FaceColor','interp','Parent',gca);
	axis([min(x) max(x)  min(y) max(y) ])	
% 	axis tight


if doline
	hold on;	
	if is12lead
		for i=2
			plot3(X(:,i),csig(:,i)./maxi+Y(1,i),1000*ones(size(X(:,i))),'k')
		end
	else
		for i=1:size(X,2)
			plot3(X(:,i),csig(:,i)./maxi+Y(1,i)+0.5,1000*ones(size(X(:,i))),'k')
		end		
	end
end
end

% caxis([-maxi*(1+1/(0.5*length(cols)-1)) maxi])
caxis([-maxi maxi])
if is12lead
	set(gca,'YTick',yt,'YTickLabel',labels,'YTickLabelMode','manual');
else
	set(gca,'YTick',y,'YTickLabel',labels,'YTickLabelMode','manual');
end

if reverseax
	set(gca,'xdir','reverse')
end

if length(yt) >1
	colorbar('peer',gca)
end
