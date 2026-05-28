function showInvCase(varargin)

global dodeprep

dodeprep=1;
bsmfile=[];
beat=[];
beat='01';
alpha=2;
layfile='nim64.mla';
type='ventricle';
if length(varargin) < 1
	error('This routine needs at least two parameters');
else
	subject=varargin{1};
	pp=2;
	while pp<=nargin
		if ischar(varargin{pp})
			key=lower(varargin{pp});
			switch key
				case 'dodeprep'
					dodeprep=varargin{pp+1};pp=pp+2;
				case 'bsmfile'
					bsmfile=varargin{pp+1};pp=pp+2;
				case 'beat'
					beat=varargin{pp+1};pp=pp+2;
				case 'layfile'
					layfile=varargin{pp+1};pp=pp+2;
				case 'mycase'
					mycase=varargin{pp+1};pp=pp+2;
				case 'type'
					type=varargin{pp+1};pp=pp+2;
				case 'alpha'
					alpha=varargin{pp+1};pp=pp+2;					
				otherwise
					error('unknown parameter');
			end
		end
	end
end

if isempty(beat)
	error('no beat specified');
	end
if isempty(bsmfile)
	bsmfile=[subject '\' subject beat '.mes'];
end

figure(1);	set(gcf,'position',[54   300   1200   600])

dirout='results\';
dirfig=['figs\' subject '\'];

if ~exist(dirfig,'dir')
	mkdir(dirfig);
end

lpass=5;
% 0.028519    0.019312


GEOM=invInit('subject',subject,'layfile',layfile,'bsm',bsmfile,'type',type,'anisotropyRatio',alpha);
GEOM=prepare_geom(GEOM,[GEOM.subject beat '_qrst.spe']);
if ~isempty(subject)
	if dodeprep
		if strfind(GEOM.type,'atria')
			mycase=[pwd '\' dirout GEOM.subject '_' beat '_alpha' num2str(alpha) GEOM.type '_deprep.mat'];
		else
			mycase=[dirout '\' GEOM.subject '_alpha2_deprep.mat'];
		end
	else
		mycase=[pwd '\' dirout GEOM.subject beat '_alpha2.mat'];
	end
end
if ~exist(mycase,'file')
	error([mycase ' does not exist']);
end
load(mycase)

REPRANGE=range(meas.repfinal);
APDRANGE=range(meas.repfinal-meas.depfinal);


% onsetqrs  =GEOM.specs(2);
% endqrs    =GEOM.specs(3);
% qrsduration=endqrs-onsetqrs+1;
if strfind(type,'atria')
	PSIREF=baselinecor(GEOM.BSM,GEOM.specs(2),GEOM.specs(3));
	PSIREF=PSIREF(:,GEOM.specs(2):GEOM.specs(3)+20);
else
	PSIREF=baselinecor(GEOM.BSM,GEOM.specs(2),GEOM.specs(5));
	PSIREF=PSIREF(:,GEOM.specs(2):GEOM.specs(5)+20);
end
t=0:size(PSIREF,2)-1;
T=ones(length(GEOM.VER),1)*t;
if strfind(type,'atria')
	Sfin=getSmode(T,meas.depfinal,meas.repfinal,GEOM.pS,[],1);	
else
	Sfin=getSmode(T,meas.depfinal,meas.repfinal,GEOM.pS,[],4);	
end
PHIfin=lowpassma(GEOM.AMA*Sfin,10);	

%%	

set(gcf,'position',[54   300   1200   600])

showAtria(GEOM.VER,GEOM.ITRI,meas.dep,'nodes',meas.foci);
saveas(gcf,[dirfig 'depinit.png']);	
	
showAtria(GEOM.VER,GEOM.ITRI,meas.depfinal,'nodes',meas.foci);
saveas(gcf,[dirfig 'depfinal.png']);	

if ~strfind(type,'atria')
	showAtria(GEOM.VER,GEOM.ITRI,meas.repfinal,'max',REPRANGE);
	saveas(gcf,[dirfig 'repfinal.png']);	

	showAtria(GEOM.VER,GEOM.ITRI,meas.rep,'max',REPRANGE);
	saveas(gcf,[dirfig 'repinit.png']);	

	showAtria(GEOM.VER,GEOM.ITRI,meas.repfinal-meas.depfinal,'max',APDRANGE);
	saveas(gcf,[dirfig 'apdfinal.png']);	

	showAtria(GEOM.VER,GEOM.ITRI,meas.rep-meas.dep,'max',APDRANGE);
	saveas(gcf,[dirfig 'apdinit.png']);	
end

clf
sigplot_p(PSIREF,'lay',GEOM.LAY,'color','b'); hold on
sigplot_p(PHIfin,'lay',GEOM.LAY,'color','r'); hold on	
saveas(gcf,[dirfig 'BSM.png']);
clf; 
if strfind(type,'atria')
	leadv16(PSIREF,PHIfin,'max',[-.2 .201],'leadsys','nim','Amplification',.1); hold on
else
	leadv16(PSIREF,PHIfin,'max',0.5,'leadsys','nim'); hold on
end
saveas(gcf,[dirfig 'ECG.png']);
