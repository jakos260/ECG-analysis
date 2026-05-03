function BipolarECG(varargin) %(t,PHI)

% plot standard leads from the nijmegen leads system


VLab=['VR';'VL';'VF';'V1';'V2';'V3';'V4';'V5';'V6'];
cols=['b','r','k','g','c','m'];
extr=1;

do9=0;
tmax=0;
marks=0;
info=[];
maxphi=0;
sampfreq=500;
delta=0;
deltamaxphi=1;
nsigs=0;
leadText=[];
markers=[];
if length(varargin) < 1
	error('This routine needs at least two parameters');
else
	pp=1;
	while pp<=nargin
		if ischar(varargin{pp})
			key=lower(varargin{pp});
			switch key
				case 'marks'
					marks=varargin{pp+1};pp=pp+2;
				case 'info'
					info=varargin{pp+1};pp=pp+2;
				case 'max'
					maxphi=varargin{pp+1};pp=pp+2;
				case 'dmax'
					deltamaxphi=varargin{pp+1};pp=pp+2;
				case 'sampfreq'
					sampfreq=varargin{pp+1};pp=pp+2;
				case 'leadtext'
					leadText=varargin{pp+1};pp=pp+2;
				case 'delta'
					delta=varargin{pp+1};pp=pp+2;
				case 'markers'
					markers=varargin{pp+1};pp=pp+2;
				otherwise
					error('unknown parameter');
			end
		else
			eval(['PHI_' num2str(pp) '=varargin{' num2str(pp) '};']);
			eval(['if size(PHI_' num2str(pp) ',2)==12,PHI_' num2str(pp) '=PHI_' num2str(pp) '''' ';end;'])
			eval(['PHI_' num2str(pp) '=cast(PHI_' num2str(pp) ',' '''' 'double' '''' ');'])
			eval(['nsigs=max(nsigs,size(PHI_' num2str(pp) ',2));']);
			if nsigs>1
				eval(['tmax=max(size(PHI_' num2str(pp) ',1)-1,tmax);']);
			else
				eval(['tmax=max(length(PHI_' num2str(pp) ')-1,tmax);']);
			end
			nmap=pp;
			pp=pp+1;
		end
	end
end
linestyl='-';
mark='none';
tmax=tmax*1000/sampfreq;
t=0:1000/sampfreq:tmax;
if delta==0
	delta=10*1000/sampfreq;
end
%% wct reference and determine maxphi

for i=1:nmap
	if maxphi==0
		eval(['maxphi=max(max(max(maxphi,(abs(PHI_' num2str(i) ')))));']);
	end
end

if ~isempty(info)
	if marks>0
		eval(['PHI_m=PHI_' num2str(marks) ';']);	
	end
end
% clf
pos=get(gcf,'Position');
maxsize=get(0,'ScreenSize');


if tmax >10000
	mark='.';
	linestyl='none';
	set(gcf,'Position',[pos(1)  min(pos(2),520) 520   520]);
else
	set(gcf,'Position',[pos(1)  min(pos(2),520) 200   520]);	
end

dm=maxphi*0.02;
linew=1.5;
fs=12;
fweight='demi';
clf

	ny=nsigs;
	for k=1:nsigs
		y=nsigs-icyc(k,nsigs);
		axes('Position',[0.05,y/ny,0.9,1/ny]);

		for i=1:nmap, 
			eval(['PHI=PHI_' num2str(i) ';']); hold on
			plot(t, PHI(:,k),cols(i),'Linewidth',linew,'linestyle',linestyl,'marker',mark);
		end
		% make sure the first signal is always seen
		eval(['PHI=PHI_' num2str(1) ';']);	hold on; 
		plot(t, PHI(:,k),[cols(1) '--'],'Linewidth',linew,'linestyle',linestyl,'marker',mark);
	
		line([0 tmax],[0 0],'Color','k','Linewidth',1);
		for i=0:delta:tmax; line([i i],[dm -dm],'Color','k','Linewidth',1); end
		axis([0,tmax,-maxphi,maxphi]); axis off
		if ~isempty(leadText)
			text(tmax*0.05,0.8*maxphi,leadText(k,:),'FontWeight',fweight,'Fontname','Verdana','FontSize',fs);
		end
		for i=1:length(markers)
			line([markers(i) markers(i)],[-maxphi+10*dm -maxphi],'Color','y','Linewidth',1);
		end
		
		

	end
	AxisLegend(t,[-maxphi maxphi],delta*2,deltamaxphi);
