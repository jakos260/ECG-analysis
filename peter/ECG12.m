function ECG12(varargin) %(t,PHI)

% plot standard leads from the nijmegen leads system


VLab=['aVR';'aVL';'aVF';'V1 ';'V2 ';'V3 ';'V4 ';'V5 ';'V6 '];
cols=['b';'r';'k';'g';'c';'m';'y'];
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
leadText='';
markers=[];
Position=[];
keep=0;
if length(varargin) < 1
	error('This routine needs at least two parameters');
else
	pp=1;
	while pp<=nargin
		if ischar(varargin{pp})
			key=lower(varargin{pp});
			switch key
				case '9leads'
					do9=varargin{pp+1};pp=pp+2;
				case 'marks'
					marks=varargin{pp+1};pp=pp+2;
				case 'do9'
					do9=varargin{pp+1};pp=pp+2;
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
				case 'keep'
					keep=varargin{pp+1};pp=pp+2;

				case 'color'
					cols=varargin{pp+1};pp=pp+2;
				case 'markers'
					markers=varargin{pp+1};pp=pp+2;
				case 'position'					
					Position=varargin{pp+1};pp=pp+2;
				otherwise
					error('unknown parameter');
			end
		else
			eval(['PHI_' num2str(pp) '=varargin{' num2str(pp) '};']);
			A=[];
			eval(['if iscell(PHI_' num2str(pp) '),A=PHI_' num2str(pp) ';end;']);
			if ~isempty(A)
				for k=1:length(A)
					eval(['PHI_' num2str(k) '=cell2mat(A(k));']);
					eval(['if size(PHI_' num2str(k) ',2)==12,PHI_' num2str(k) '=PHI_' num2str(k) '''' ';end;'])
					eval(['PHI_' num2str(k) '=cast(PHI_' num2str(k) ',' '''' 'double' '''' ');'])
					eval(['nsigs=max(nsigs,size(PHI_' num2str(k) ',1));']);
					eval(['tmax=max(size(PHI_' num2str(k) ',2)-1,tmax);']);
				end
				nmap=length(A);
				linew=1;
				cols='b';
				pp=pp+1;
			else
				eval(['if size(PHI_' num2str(pp) ',2)==12,PHI_' num2str(pp) '=PHI_' num2str(pp) '''' ';end;'])
				eval(['PHI_' num2str(pp) '=cast(PHI_' num2str(pp) ',' '''' 'double' '''' ');'])
				eval(['nsigs=max(nsigs,size(PHI_' num2str(pp) ',1));']);
				eval(['tmax=max(size(PHI_' num2str(pp) ',2)-1,tmax);']);
				nmap=pp;
				pp=pp+1;
			end
		end
	end
end
linestyl='-';
mark='none';
tmax=tmax*1000/sampfreq;
t=0:1000/sampfreq:tmax;
if delta==0
	delta=100*length(t)/sampfreq;
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
	else
		marks=1;			
	end
elseif marks>0
	eval(['PHI_m=PHI_' num2str(marks) ';']);
	info=getAmplInfos(PHI_m);
end

% clf
pos=get(gcf,'Position');
maxsize=get(0,'ScreenSize');
if ~isempty(Position)
	set(gcf,'Position',Position);
elseif do9
	set(gcf,'Position',[pos(1)  min(pos(2),520) 1060*2/3   300]);
else
	set(gcf,'Position',[pos(1)  min(pos(2),maxsize(4)-1.33333*520-150) 1060   700]);
end
if tmax >10000
	mark='.';
	linestyl='none';
	if do9
		set(gcf,'Position',[pos(1)  min(pos(2),520) 1400   520]);
	else
		set(gcf,'Position',[pos(1)  min(pos(2),maxsize(4)-1.33333*520-150) 1400   520]);
	end
end

dm=maxphi*0.02;
linew=1.5;
fs=12;
fweight='demi';
if ~keep,clf;end


if do9
	ny=3;nx=3;
	x=0;
else
	x=0;
	ny=3;nx=4;
	textI=['I  ';'II ';'III'];
	for k=1:3
		y=3-icyc(k,3);
		axes('Position',[x/nx+0.025,y/ny,1/nx-0.05,1/ny]);
		for i=1:nmap, 
			eval(['PHI=PHI_' num2str(i) ';']); hold on
			if ~isempty(PHI)
				if k==1
					plot(t, PHI(1,:),cols(icyc(i,size(cols,1)),:),'Linewidth',linew,'linestyle',linestyl,'marker',mark);
				elseif k==2
					plot(t, PHI(2,:),cols(icyc(i,size(cols,1)),:),'Linewidth',linew,'linestyle',linestyl,'marker',mark);
				else
					plot(t, PHI(3,:),cols(icyc(i,size(cols,1)),:),'Linewidth',linew,'linestyle',linestyl,'marker',mark);				
				end
			end
		end
		
		axis([0,tmax,-maxphi,maxphi]); axis off
		if ~keep
			% make sure the first signal is always seen
			eval(['PHI=PHI_' num2str(1) ';']);	hold on; 
			if ~isempty(PHI)
				if k==1
					plot(t, PHI(1,:),[cols(1,1) '--'],'Linewidth',linew,'linestyle',linestyl,'marker',mark);
				elseif k==2
					plot(t, PHI(2,:),[cols(1,1) '--'],'Linewidth',linew,'linestyle',linestyl,'marker',mark);
				else
					plot(t, PHI(3,:),[cols(1,1) '--'],'Linewidth',linew,'linestyle',linestyl,'marker',mark);				
				end
			end			
			line([0 tmax],[0 0],'Color','k','Linewidth',1);
			for i=0:delta:tmax; line([i i],[dm -dm],'Color','k','Linewidth',1); end

			text(tmax*0.05,0.8*maxphi,textI(k,:),'FontWeight',fweight,'Fontname','Verdana','FontSize',fs);	
			for i=1:length(markers)
				line([markers(i) markers(i)],[-10*dm 10*dm ],'Color','y','Linewidth',1);
			end
		end
	end
	x=1;
end

for k=1:9
	y=3-icyc(k,3);
	axes('Position',[x/nx+0.025,y/ny,1/nx-0.05,1/ny]);
	if mod(k,ny)==0, x=x+1; end
	for i=1:nmap, 
		eval(['PHI=PHI_' num2str(i) ';']); hold on
		if ~isempty(PHI)
			plot(t, PHI(k+3,:),cols(icyc(i,size(cols,1)),:),'Linewidth',linew,'linestyle',linestyl,'marker',mark);
		end
	end
	
	axis([0,tmax,-maxphi,maxphi]); axis off
	if ~keep
	% make sure the first signal is always seen
		eval(['PHI=PHI_' num2str(1) ';']);	hold on; 
		if ~isempty(PHI)
			plot(t, PHI(k+3,:),[cols(1,1) '--'],'Linewidth',linew,'linestyle',linestyl,'marker',mark);
		end			
		line([0 tmax],[0 0],'Color','k','Linewidth',1);
		for i=0:delta:tmax; line([i i],[dm -dm],'Color','k','Linewidth',1); end
		axis([0,tmax,-maxphi,maxphi]); 
		text(tmax*0.05,0.8*maxphi,VLab(k,:),'FontWeight',fweight,'Fontname','Verdana','FontSize',fs);
		for i=1:length(markers)
			line([markers(i) markers(i)],[-10*dm 10*dm ],'Color','y','Linewidth',1);
		end
	end
end
if ~keep
	AxisLegend(t,[-maxphi maxphi],delta*2,deltamaxphi);
end
	axis([0,tmax,-maxphi,maxphi]); axis off