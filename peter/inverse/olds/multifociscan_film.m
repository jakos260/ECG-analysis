function [bestfoci,bestdep,notchPot,outp]=multifociscan_film(GEOM,leads,scanmode)

% use all foci found each round

% date:30-04-2009
% identification of one or more focal points of depolarization,
% Peter van Dam
% global lpass
% SCAN SPECIFIC PARAMETERS
% the parameters shown below are the only remaining ones 

usecor=1;
global dodeprep
global FF
global GG
%%
if scanmode>1
	if dodeprep
		notchPot=sqrAmpl(GEOM,GEOM.specs(3)+10:GEOM.specs(3)+40,leads);
	else
		notchPot=1-sqrAmpl(GEOM,GEOM.specs(3)+10:GEOM.specs(3)+40,leads);
	end
else
	notchPot=zeros(size(GEOM.DIST2W(:,1)));
end


%% Scan for foci
qrsduration=GEOM.specs(3)-GEOM.specs(2)+1;
dep=[];bestfoci=[];outp=[];bestdep=dep;bestcor=-100;bestrd=100;
startTime=clock;
k=0;
while  length(unique(bestfoci))<100 
	k=k+1;
	tic
	[nofocus,dep,foci,bestcor,bestrd,bestshift]=fociscan(GEOM,usecor,dep,leads,scanmode,notchPot,bestcor);	
	dt=ceil(toc);
	if nofocus, 
		break;
	end
	if usecor &&~isempty(outp)&& round(10000*bestcor) <= round(10000*max(outp(:,1))),	break;	end	 
	if ~usecor&&~isempty(outp)&& bestrd >= max(outp(end,2)),break;	end	 

	bestdep=dep;bestfoci=[bestfoci foci];
	outp=[outp;[bestcor bestrd bestshift max(dep)]];
	nrClust=nrClusters(GEOM,unique(bestfoci));
	disp(['nr foci/clusters: ' num2str(length(unique(bestfoci)) ,3) '  ' num2str(nrClust,3) ...
		  '   QRS duration /sim: ' num2str(qrsduration,3) '  ' num2str(max(dep),3)...
		  '   corr/rd/init delay: ' num2str(bestcor,3) ' ' num2str(bestrd,3) ' ' num2str(bestshift,2)...
		  '  ' datestr(datenum(clock)-datenum(startTime),'HH,MM.SS')])
	if k==2
	hdep=figure(13);set(gcf,'color','w');showAtria(GEOM.VER,GEOM.ITRI,bestdep,'nodes',bestfoci,'view',[-5 50],'endo',GEOM.endoVER);drawnow		
	else
	hdep=figure(13);set(gcf,'color','w');showAtria(GEOM.VER,GEOM.ITRI,bestdep,'nodes',bestfoci,'view',[-55 50],'endo',GEOM.endoVER);drawnow				
	end
% 	annotation('textbox',[0. 0.9 0.3 0.05],'string',['# ' num2str(k) '       Correlation: ' num2str(bestcor,2) '   rd:' num2str(bestrd,2) ],'edgecolor','none','FitHeightToText','on','Fontsize',24');
	annotation('textbox',[0. 0.9 0.3 0.05],'string',['# ' num2str(k) ],'edgecolor','none','FitHeightToText','off','Fontsize',12');	
	saveas(hdep,['figs\present\multifocal_base' num2str(k) '.png']);
	
	if ~isempty(FF)
		for pp=1:min(dt,1)
			FF=[FF getFrame(hdep)];
		end
	else
		FF=getFrame(hdep);
	end
	if k==1
		dep1=bestdep;
	elseif k==2
	hdep=figure(13);set(gcf,'color','w');showAtria(GEOM.VER,GEOM.ITRI,dep1,'nodes',bestfoci,'view',[-5 50],'endo',GEOM.endoVER);drawnow		
	annotation('textbox',[0. 0.9 0.3 0.05],'string','# 1' ,'edgecolor','none','FitHeightToText','off','Fontsize',12');	
	saveas(hdep,['figs\present\multifocal_base_extra1'  '.png']);
	end	
	qrsduration=GEOM.specs(3)-GEOM.specs(2)+1;
	PSIREF=GEOM.BSM(:,GEOM.specs(2):GEOM.specs(5));
	PSIREF=baselinecor(zeromean(PSIREF));
	PSIREF=PSIREF(leads,1:qrsduration);
	t=0:size(PSIREF,2)-1;
	T=ones(length(GEOM.VER),1)*t;
	PSIA=baselinecor(lowpassma(GEOM.AMA*getSmode(T,bestdep,[],GEOM.pS,[],scanmode),5));
	hecg=figure(114);set(gcf,'color','w','PaperPositionMode','auto');
	leadv16(PSIREF,PSIA , 'leadsys',GEOM.leadname,'paperspeed',200,'max',[-3 3]	);drawnow			
% 	annotation('textbox',[0. 0.95 0.01 0.05],'string',['# ' num2str(iter) '       rd:' num2str(TST.rd,2) ],'edgecolor','none','FitHeightToText','on','Fontsize',14');		
	annotation('textbox',[0. 0.91 0.12 0.1],'string',sprintf('%s \n%s \n%s',['# ' num2str(k)],['R:' num2str(bestcor,2)],['rd:' num2str(bestrd,2) ]),'edgecolor','none','FitHeightToText','off','Fontsize',12,'fontweight','normal','BackgroundColor','w','color','k');
	saveas(hecg,['figs\present\multifocal_ECG' num2str(k) '.png']);

	if ~isempty(GG)
		for pp=1:min(dt,1)
			GG=[GG getFrame(gcf)];
		end
	else
		GG=getFrame(gcf);
	end
end
% [bestdep,outp(end,1),outp(end,2)]=delayscan(GEOM,bestdep,leads,scanmode,notchPot);	
% 
% disp('final')
% disp(['nr foci/clusters: ' num2str(length(unique(bestfoci)) ,3) '  ' num2str(nrClust,3) ...
% 	  '   QRS duration /sim: ' num2str(qrsduration,3) '  ' num2str(max(dep),3)...
% 	  '   corr/rd/init delay: ' num2str(outp(end,1),3) ' ' num2str(outp(end,2),3) ' ' num2str(bestshift,2)...
% 	  '  ' datestr(datenum(clock)-datenum(startTime),'HH,MM.SS')])

%% =======================================================================
function [nofocus,bestdep,foci,bestcor,bestrd,bestshift]=fociscan(GEOM,usecor,initdep,leads,scanmode,notchPot,initcor)

global lpass
% fixed paramters

MAX_PURK_SHIFT=0.4;
MAX_MYO_SHIFT=0.9;
MAX_MYO_VS=0.8;
if strfind(GEOM.type,'atria')
	MAX_MYO_VS=0.9;
	MAX_MYO_SHIFT=0.8;
end

%% prepare 
% ECG signals are only used between start QRS and end Twave
qrsduration=GEOM.specs(3)-GEOM.specs(2)+1;
PSIREF=GEOM.BSM(:,GEOM.specs(2):GEOM.specs(5));
PSIREF=baselinecor(zeromean(PSIREF));
PSIREF=PSIREF(leads,1:qrsduration);
normphi=norm(PSIREF,'fro');

rep=100*ones(size(GEOM.VER,1),1)+GEOM.specs(4);%only needed if scanmode ~=1
t=0:size(PSIREF,2)-1;
T=ones(length(GEOM.VER),1)*t;
AMA=zeromean(GEOM.AMA(leads,:));

% init
foci=0;	
nofocus=0;bestdep=initdep;bestcor=-1;bestrd=10;bestshift=-1;
cors=-1*ones(size(GEOM.VER,1),1);rds=10*ones(size(GEOM.VER,1),1);
deps=zeros(size(GEOM.ADJ));		% all depolarization sequences
VS=zeros(length(GEOM.ADJ),1);
% shifts
if isempty(initdep)
	shift=0*ones(size(GEOM.DIST2W(:,1)));
else
	% initial shifts are MAX_PURK_SHIFT% of the previous depolarization time 
	% for purkinje nodes and MAX_MYO_SHIFT % for myocardial nodes. In case
	% the previous depolarization duration is longer than the QRS complex 
	% also the purkinje nodeshift is allowed on myocardial nodes. Most
	% probably there is more than one focus in this situation. The reason
	% for making a diffrence in myocardial nodes and purkinje nodes is that
	% under almost all circumcantces the activation wave ends on myocardial
	% nodes. Optimizing the depolarization sequence without changing this
	% latest site of activation causes only local small minima, resulting
	% in an alternative way of inverse timing estimation, without teh
	% posibility to delay activation!
	if abs(max(initdep)-qrsduration)<1e-10 %ms
		shift(GEOM.purkinjever==1)=initdep(GEOM.purkinjever==1)*MAX_PURK_SHIFT;
		shift(GEOM.purkinjever==0)=initdep(GEOM.purkinjever==0)*MAX_MYO_SHIFT;
	else
		shift=initdep*MAX_PURK_SHIFT;
	end
	shift(shift<1)=0; % prevent suboptimalization
end
%% scan for foci -----------------------------------------------
for inode=1:length(GEOM.VER)
	% these nodes are already focus and cannot be activated earlier
	if ~isempty(initdep)&& shift(inode)<1, continue;end 
	dep=shift(inode)+GEOM.DIST2W(:,inode);
	if ~isempty(initdep) % no previous focus determined
		if ~GEOM.purkinjever(inode) %myocarial node
			% determine vs (estimation of the longitudinal velocity) for a
			% myocardial node limiting the velocity. 
			vs=min(max(min(dep,initdep))/qrsduration,MAX_MYO_VS);	
			dep=min(dep/vs,initdep);
		else % purkinje node
			dep=min(dep,initdep);				
		end
		vs=max(dep)/(qrsduration);		
		dep=dep/vs;	
		if sum(dep==initdep)/length(dep)>0.95, continue;end;% nothing to do because not enough nodes are affected
	else % first run
		vs=max(dep)/(qrsduration);
		% if inode is a myocardial node the maximum velocity overrules the 
		% fact that the qrs duration and max(dep) should be equal
		if ~GEOM.purkinjever(inode),vs=min(vs,MAX_MYO_VS);	end
		dep=dep/vs;			
		VS(inode)=vs;
	end
	deps(:,inode)=dep;
	PSIA=lowpassma(baselinecor(AMA*getSmode(T,dep,rep,GEOM.pS,notchPot,scanmode)),lpass);
% 	PSIA=(lowpassma(AMA*getSmode(T,dep,rep,GEOM.pS,notchPot,scanmode),lpass));	
	COR=corrcoef(PSIA,PSIREF);if COR>cors(inode),cors(inode)=COR(2,1);end	
	RD=norm(PSIREF-PSIA,'fro')/normphi;if RD<rds(inode),rds(inode)=RD;end
end
%%
if usecor

	if isempty(initdep) && qrsduration<=110
		% if the qrs duration is relatively short teh first focus should be
		% found in the left ventricle
		cors(GEOM.Rfreewallver==1)=-1;
	end
	% mcors is experimental. 1 vogel (focus) maakt nog geen lente)	
	mcors=cors;
	for i=1:length(cors)
		a=find(GEOM.ADJsurf(i,:)<20 &GEOM.ADJsurf(i,:)>0);
		a(cors(a)<cors(i)-0.1)=[];
		mcors(i)=mean(cors([i a]));
	end
	id=5; % mcors
	if max(mcors)<initcor && max(cors)>initcor,
		id=2;
	end
	id=2;
	% detemine the node with the highest correlation
	A=[(1:length(cors))' cors rds VS mcors];A=sortrows(A,id);A=A(end:-1:1,:);
	if abs(A(1,2)-mean(A(1:3,2)))<1e-10,return;	end
	A(A(:,2)<=0,:)=[];	
	if isempty(A)
		return;
	else
		A=A(1,:);
		bestdep=deps(:,A(1));
		foci=A(1);
		bestcor=A(id);
		bestrd=A(3);
		bestshift=shift(A(1));
	end
else
	A=[(1:length(cors))' cors rds];A=sortrows(A,3);
	A(A(:,3)>=10,:)=[];
	if isempty(A)
		nofocus=1;
		return;
	else
		A=A(1,:);
		bestdep=deps(:,A(1));
		foci=A(1);
		bestcor=A(2);
		bestrd=A(3);
		bestshift=shift(A(1));
	end
end      

%% delay scan-----------------------------------------------
% experiment, probably not usefull. 
function [bestdep,bestcor,bestrd]=delayscan(GEOM,initdep,leads,scanmode,notchPot)

global lpass;

bestdep=initdep;
qrsduration=GEOM.specs(3)-GEOM.specs(2)+1;
PSIREF=GEOM.BSM(:,GEOM.specs(2):GEOM.specs(5));
PSIREF=baselinecor(zeromean(PSIREF));
PSIREF=PSIREF(leads,1:qrsduration);
normphi=norm(PSIREF,'fro');

rep=100*ones(size(GEOM.VER,1),1)+GEOM.specs(4);%only needed if scanmode ~=1
t=0:size(PSIREF,2)-1;
T=ones(length(GEOM.VER),1)*t;
AMA=zeromean(GEOM.AMA(leads,:));

PSIA=baselinecor(lowpassma(AMA*getSmode(T,initdep,rep,GEOM.pS,notchPot,scanmode),lpass));
COR=corrcoef(PSIA,PSIREF);
bcor=COR(2,1);
bestrd=norm(PSIREF-PSIA,'fro')/normphi;
%%
keepdep=initdep;
DELAYFACT=0.005;
keepdelay=zeros(size(keepdep));
for i=1:10
	for inode=1:length(GEOM.VER)
		dep=keepdep+keepdelay;
		delay=keepdelay;
		nodes=[inode find(GEOM.DISTsurf(inode,:)<30)];
		delay(nodes)=((-1)^i)*dep(nodes)*DELAYFACT;
		dep=dep+delay;
		dep(dep>qrsduration)=qrsduration;
		dep(dep<0)=0;	
		dep(nodes)=max(dep(nodes),0.1*keepdep(nodes));
		dep(nodes)=min(dep(nodes),1.1*keepdep(nodes));
		PSIA=baselinecor(lowpassma(AMA*getSmode(T,dep,rep,GEOM.pS,notchPot,scanmode),lpass));
		COR=corrcoef(PSIA,PSIREF);COR=COR(2,1);	RD=norm(PSIREF-PSIA,'fro')/normphi;
		RD=norm(PSIREF-PSIA,'fro')/normphi;
		if RD<bestrd && round(100*COR)/100>=round(100*bcor)/100
			bestrd=RD;bcor=COR;	bestdep=dep;
			keepdelay(nodes)=dep(nodes)-keepdep(nodes);		
		end
	end
end
bestcor=bcor;
% figure(14);showAtria(GEOM.VER,GEOM.ITRI,(keepdep-bestdep)./keepdep,'max',
% [-0.1 0.1]);