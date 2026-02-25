function feature12lead(varargin)

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
	dodiff=0;
	dozmap=0;
	xlabels=[];
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
				case 'markers'
					markers=varargin{pp+1};pp=pp+2;
				case 'episodes'
					episodes=varargin{pp+1};pp=pp+2;
				case 'dodiff'
					dodiff=varargin{pp+1};pp=pp+2;
				case 'dozmap'
					dozmap=varargin{pp+1};pp=pp+2;
				case 'xlabel'
					xlabels=varargin{pp+1};pp=pp+2;
					
				otherwise
					error('unknown parameter');
			end
		end
	end
end

if ~isempty(episodes)
	type='episodes';
	markers=episodes;
else
	type='markers';
end
if level<0
	maxtype='max';
else
	maxtype='level';
	maxi=level;
end
% labs=['I  ';'II ';'III';'aVR';'aVL';'aVF';'V1 ';'V2 ';'V3 ';'V4 ';'V5 ';'V6 '];
% ind=1:12;
% labs=['I  ';'aVF';'aVR';'V1 ';'V2 ';'V3 ';'V4 ';'V5 ';'V6 ';'aVL';'III';'II '];
% ind=[1 5 4 7 8 9 10 11 12 6  3 2];    
labs=['aVR';'aVL';'I  ';'V1 ';'V2 ';'V3 ';'V4 ';'V5 ';'V6 ';'II ';'aVF';'III'];
ind=[4 5 1 7 8 9 10 11 12 2 6 3];    

clf
for k=1:12
	axes('Position',[0.05 k/13 0.8 1/15])
	featuremap(sig(:,ind(13-k)),maxtype,maxi,type,markers,'labels',labs(13-k,:),'dodiff',dodiff,'dozmap',dozmap);
	box on
	if k>1
		set(gca,'XTick',[]);
	else
		if ~isempty(xlabels)
			set(gca,'XTick',markers,'XTicklabel',xlabels);
		end
		axis on;
	end
end

if dozmap
	maxi=3;
	cols=loadmat('sigmamap.mcm');
else
	cols=loadmat('mymap_basic.mcm');
end
colormap(cols)



axes('Position',[0.88 0.05 0.09 0.9])
colorbar('East')
axis off
% caxis([-maxi*(1+1/(0.5*length(cols)-1)) maxi])
caxis([-maxi maxi])

