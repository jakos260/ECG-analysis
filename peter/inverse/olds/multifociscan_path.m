function [bestfoci,bestdep,outp]=multifociscan_path(GEOM,anisotropyRatio,leads)

% use all foci found each round

% date:15032009
% identification of one or more focal points of depolarization,
% Peter van Dam
% global lpass
% SCAN SPECIFIC PARAMETERS
% the parameters shown below are the only remaining ones 

scanmode=1;
usecor=1;


%% anisotropic matrix calculation
if anisotropyRatio==1
	GEOM.ADJ2W=GEOM.ADJ;
	DIST=GEOM.DIST;
else
	GEOM.ADJ2W=calcAnisADJ(GEOM,anisotropyRatio);	
	if strcmp(GEOM.type,'_atria')
		if exist(['./temp/' GEOM.subject num2str(anisotropyRatio) GEOM.type '.dist'],'file')
			DIST=loadmat(['./temp/' GEOM.subject num2str(anisotropyRatio) GEOM.type '.dist']);
		else
			DIST=graphdist(GEOM.ADJ2W);
			savemat(['./temp/' GEOM.subject num2str(anisotropyRatio) GEOM.type '.dist'],DIST);
		end
	elseif exist(['./temp/' GEOM.subject num2str(anisotropyRatio) '.dist'],'file')
		DIST=loadmat(['./temp/' GEOM.subject num2str(anisotropyRatio) '.dist']);
	else
		DIST=graphdist(GEOM.ADJ2W);
		savemat(['./temp/' GEOM.subject num2str(anisotropyRatio) '.dist'],DIST);
	end
end

%% Scan for foci
qrsduration=GEOM.specs(3)-GEOM.specs(2)+1;
dep=[];bestfoci=[];outp=[];
startTime=clock;

while  length(unique(bestfoci))<100 
	[nofocus,dep,foci,bestcor,bestrd,bestshift]=fociscan(GEOM,usecor,dep,DIST,leads,scanmode,bestfoci);	
	if nofocus, 
		break;
	end
	if usecor &&~isempty(outp)&& bestcor <= max(outp(:,1)),	break;	end	 
	if ~usecor&&~isempty(outp)&& bestrd >= max(outp(end,2)),break;	end	 	
	outp=[outp;[bestcor bestrd bestshift max(dep)]];
	bestdep=dep;bestfoci=[bestfoci foci];
	nrClust=nrClusters(GEOM,unique(bestfoci));
	disp(['nr foci/clusters: ' num2str(length(unique(bestfoci)) ,3) '  ' num2str(nrClust,3) ...
		  '   QRS duration /sim: ' num2str(qrsduration,3) '  ' num2str(max(dep),3)...
		  '   corr/rd/init delay: ' num2str(bestcor,3) ' ' num2str(bestrd,3) ' ' num2str(bestshift,2)...
		  '  ' datestr(datenum(clock)-datenum(startTime),'HH,MM.SS')])
	figure(13);showAtria(GEOM.VER,GEOM.ITRI,bestdep,'nodes',bestfoci);drawnow	
% 	if strfind(GEOM.type,'atria'), usecor=0;end
end

%% =======================================================================
function [nofocus,bestdep,foci,bestcor,bestrd,bestshift]=fociscan(GEOM,usecor,initdep,DIST2W,leads,scanmode,initfoci)

global lpass
% fixed paramters
MAX_MYO_VS=0.8;
DELTA_SHIFT=0.8;
MAX_MYO_SHIFT=0.8;
if strfind(GEOM.type,'atria')
	MAX_MYO_VS=.9;
	MAX_MYO_SHIFT=0.7;
end

%% prepare 
% ECG signals are only used between start QRS and end Twave
qrsduration=GEOM.specs(3)-GEOM.specs(2)+1;
PSIREF=GEOM.BSM(:,GEOM.specs(2):GEOM.specs(5));
PSIREF=baselinecor(PSIREF);
PSIREF=PSIREF(leads,1:qrsduration);
normphi=norm(PSIREF,'fro');

rep=100*ones(size(GEOM.VER,1),1)+GEOM.specs(4);%only needed if scanmode ~=1
t=0:qrsduration-1;
T=ones(length(GEOM.VER),1)*t;
AMA=GEOM.AMA(leads,:);

% init
nofocus=0;bestdep=initdep;bestcor=-1;bestrd=10;bestshift=-1;
cors=-1*ones(size(GEOM.VER,1),1);rds=10*ones(size(GEOM.VER,1),1);
deps=zeros(size(GEOM.ADJ));		% all depolarization sequences

% shifts
if isempty(initdep)
	shift=zeros(size(DIST2W(:,1)));
else
	% initial shifts are DELTA_SHIFT% of the previous depolarization time 
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
	if max(initdep)<=qrsduration
		shift(GEOM.purkinjever==1)=initdep(GEOM.purkinjever==1)*DELTA_SHIFT;
		shift(GEOM.purkinjever==0)=initdep(GEOM.purkinjever==0)*MAX_MYO_SHIFT;
	else
 		shift=initdep*DELTA_SHIFT;
	end
	shift(shift<1)=0; % prevent suboptimalization
end

%% scan for foci -----------------------------------------------
for inode=1:length(GEOM.VER)
	% these nodes are already focus an cannot be activated earlier
	if ~isempty(initdep)&& shift(inode)==0, continue;end 
	dep=DIST2W(:,inode);
	if ~isempty(initdep) 
		if ~GEOM.purkinjever(inode) %myocarial node
			% determine vs (estimation of the longitudinal velocity) for a
			% myocardial node limiting the velocity. 
% 			vs=min(max(min(dep+shift(inode),initdep))/qrsduration,MAX_MYO_VS);	
% 			vs=max(vs,0.9);
% 			
			vs=MAX_MYO_VS;
% 			fi=find(DIST2W(inode,initfoci(initdep(initfoci)<shift(inode)))==min(DIST2W(inode,initfoci(initdep(initfoci)<shift(inode)))));
% 			if strfind(GEOM.type,'atria')
% 				fi=1;
% 			else
				fi=find(DIST2W(inode,initfoci)==min(DIST2W(inode,initfoci)));
% 			end
% 			vs=min(min(DIST2W(inode,initfoci(1)))/qrsduration,MAX_MYO_VS);
			
			if GEOM.DISTsurf(inode,initfoci(fi)) < 1.3*GEOM.DIST(inode,initfoci(fi))
				froute=getroute(GEOM.ADJsurf,initfoci(fi),inode);
				deltash=shift(inode)-initdep(initfoci(fi));
				sh=GEOM.DISTsurf(froute);sh=sh./max(sh)*deltash+initdep(initfoci(fi));
			else
				froute=getroute(GEOM.ADJ,initfoci(fi),inode);
				deltash=shift(inode)-initdep(initfoci(fi));
				sh=GEOM.DIST(froute);sh=sh./max(sh)*deltash+initdep(initfoci(fi));
			end
			sh=sh(end:-1:1);
			depr=DIST2W(:,froute)/vs;
			for i=1:length(froute)
				depr(:,i)=depr(:,i)+sh(i);
			end
% 			vs=min(min(DIST2W(inode,initfoci))/qrsduration,MAX_MYO_VS);
			dep=min([depr initdep],[],2);
		else % purkinje node
			dep=min(dep+shift(inode),initdep);				
		end
		vs=max(dep)/(qrsduration);		
		dep=dep/vs;	
% 		if all(dep==initdep),continue;end % nothing to do
		if sum(dep==initdep)/length(dep)>0.99, continue;end;% nothing to do
	else % first run
		vs=max(dep)/(qrsduration);
		% if inode is a myocardial node the maximum velocity overrules the 
		% fact that the qrs duration and max(dep) should be equal
		if ~GEOM.purkinjever(inode),vs=min(vs,MAX_MYO_VS);	end
		dep=dep/vs;			
	end
	deps(:,inode)=dep;
	PSIA=lowpassma(AMA*getSvplateau(T,dep,rep,GEOM.pS,[],scanmode),lpass);
	COR=corrcoef(PSIA,PSIREF);if COR>cors(inode),cors(inode)=COR(2,1);end	
	RD=norm(PSIREF-PSIA,'fro')/normphi;if RD<rds(inode),rds(inode)=RD;end
end
%%
if usecor
	A=[(1:length(cors))' cors rds];A=sortrows(A,2);A=A(end:-1:1,:);
	A(A(:,2)<=0,:)=[];	
	if isempty(A)
		nofocus=1;
		bestdep=[];
		foci=0;
		bestcor=-100;
		bestrd=100;
		bestshift=-1;		
		return;
	else
		A=A(1,:);
		bestdep=deps(:,A(1));
		foci=A(1);
		bestcor=A(2);
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
    

