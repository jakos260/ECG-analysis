function [bestfoci,bestdep,outp,vs]=mymultifociscan(GEOM,anisotropyRatio,leads,maxScans)
% function [foci,bestdep,GEOM]=multifociscan(GEOM,vs,mode,anisotropyRatio,leads)

% use all foci found each round

% date:161208

% identification of one or more focal points of depolarization,

% SCAN SPECIFIC PARAMETERS
% the parameters shown below are the only remaining ones 
scanmode=1;
usecor=1;
onsetqrs  =GEOM.specs(2);
endqrs    =GEOM.specs(3);
qrsduration=endqrs-onsetqrs+1;

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

vs=round(100*min(max(DIST))/(qrsduration))/100;
disp(['vs: ' num2str(vs) '  vt: ' num2str(vs/anisotropyRatio)]);
DIST=DIST/vs;

%% Scan for foci
% The actual foci scan. Added foci that yield a lower correaltation than 
% the maximum correlation - 'corrange' are removed to limit the calculation
% time. 4 succesive foci scans are done to find foci.
foci=[];dep=[];bestfoci=[];
% nrClust=0;
usetime=round(qrsduration*1);
outp=[];
maxshifttime=qrsduration-75;
tbestcor=-1;
while  length(foci)<maxScans%  nrClust<maxScans 
	[nofocus,dep,foci,bestcor,bestrd,bestshift]=focalscan(GEOM,usecor,maxshifttime,usetime,dep,foci,DIST,leads,scanmode);
	if bestcor<tbestcor || nofocus || max(dep)<qrsduration,	break;	end	 
    nrClust=nrClusters(GEOM,foci);	
	disp(['nr foci/clusters: ' num2str(length(foci) ,3) '  ' num2str(nrClust,3) ...
		  '   QRS duration /sim: ' num2str(qrsduration,4) '  ' num2str(max(dep),4)...
		  '   correlation/rd/init delay: ' num2str([bestcor bestrd bestshift],3)])
	tbestcor=bestcor;
	figure(11);showAtria(GEOM.VER,GEOM.ITRI,dep,'nodes',foci);drawnow	
	outp=[outp;[bestcor bestrd bestshift max(dep)]];
	bestdep=dep;
	bestfoci=foci;
	usetime=min(qrsduration, usetime+10);
end
% [dep,bestshift]=optimizeShifts(GEOM,DIST,bestfoci,outp(:,3),leads);
% bestdep=dep;
%% =======================================================================
function [nofocus,bestdep,foci,bestcor,bestrd,bestshift]=focalscan(GEOM,usecor,maxshifttime,minusetime,initdep,initfoci,DIST2W,leads,scanmode)

nfoci=1;
nofocus=0;

global lpass
% global clusterDist

% scan foci ; 
onsetqrs  =GEOM.specs(2);
endqrs    =GEOM.specs(3);
qrsduration=endqrs-onsetqrs+1;
minusetime=min(qrsduration,minusetime);
rep=100*ones(size(GEOM.VER,1),1)+GEOM.specs(4);
% only p(1), setting the upslope of the TMP, is effective in this
% program (mode of gets=1)
t=0:qrsduration-1;
T=ones(length(GEOM.VER),1)*t;
%%
fociver=ones(size(GEOM.purkinjever));
% fociver(initfoci)=0;
%%
PSIREF=GEOM.BSM(:,GEOM.specs(2):end);
PSIREF=baselinecor(PSIREF);
PSIREF=PSIREF(:,1:qrsduration);
if scanmode==1
	PSIREF=baselinecor(PSIREF);
end

normphi=norm(PSIREF(leads,:),'fro');
AMA=zeromean(GEOM.AMA(leads,:));

%% scan for foci -----------------------------------------------
cors=-1*ones(size(GEOM.VER,1),1);shiftRD=cors;shiftCOR=cors;
rds=100*ones(size(GEOM.VER,1),1);


deps=zeros(size(GEOM.VER,1));
maxshiftT=maxshifttime;
steps=min(10,(maxshiftT)/5);
maxshift=maxshiftT;

for inode=1:length(GEOM.VER),
	if fociver(inode)
        depi=DIST2W(:,inode);
		depi(inode)=0;
		for shift=0:steps:maxshift
			dep=depi+shift;
			if ~isempty(initdep),dep=min([initdep dep],[],2);end
			if max(dep)< qrsduration,continue;end
			S=getSvplateau(T,dep,rep,GEOM.pS,[],scanmode);
			PSIA=lowpassma(AMA*S,lpass); 
			usetime=min(qrsduration,minusetime);
	  		COR=corrcoef(PSIA(:,1:usetime),PSIREF(leads,1:usetime)); COR=COR(2,1);
			RD=norm(PSIREF(leads,1:usetime)-PSIA(:,1:usetime),'fro')/normphi;
% 			if RD>rdinit, continue;end;
			if RD<rds(inode),rds(inode)=RD;shiftRD(inode)=shift;end
			if COR>cors(inode),
% 				figure(10);	showAtria(GEOM.VER,GEOM.ITRI,dep,'nodes',inode);
				cors(inode)=COR;shiftCOR(inode)=shift;
			end			
			deps(:,inode)=dep;
		end
	end
end
% cors(shiftCOR==maxshift)=-1;
% rds(shiftRD==maxshift)=10;
% if (~usecor && sum(rds>=10)>length(shiftRD)-1)
% 	usecor=1;
% end	
if (usecor && sum(cors==-1)==length(cors))
	nofocus=1;
	return;
end

%% Determine new foci based on the scans

if usecor
	cors(GEOM.purkinjever==0)=cors(GEOM.purkinjever==0)-0.1;
	A=[(1:length(cors))' cors];A=sortrows(A,2);A=A(end:-1:1,:);
	foci=A(1,1);
	cors(GEOM.purkinjever==0)=cors(GEOM.purkinjever==0)+0.1;	
% 	A(A(:,2)<=0,:)=[];
% 	[nr,id]=nrClusters(GEOM,A(:,1));
% 	C=(1:max(id))';C=[C zeros(size(C)) zeros(size(C))];
% 	for i=1:max(id)
% 	   C(i,2)=mean(A(id==i,2));
% 	  C(i,3)=length(find(id==i));    
% 	end
% 	C=sortrows(C,2);
% 	for i=size(C,1):-1:1
% 	   if C(i,3)>1,  break; end
% 	end
% 	Aselect=A(id==C(i,1),:);
% 	foci=Aselect(1:min(nfoci,size(Aselect,1)),1);
else
	B=[(1:length(rds))' rds];B=sortrows(B,2); B(B(:,2)==0,:)=[];
	[nr,id]=nrClusters(GEOM,B(:,1));
	C=(1:nr)';C=[C zeros(size(C)) zeros(size(C))];
	for i=1:nr
		C(i,2)=mean(B(id==i,2));
	    C(i,3)=length(find(id==i));    
	end
	C=sortrows(C,2);
	for i=1:length(C)
		if C(i,3)>1, break;   end
	end
	Bselect=B(id==C(i,1),:);
	foci=Bselect(1:min(nfoci,length(Bselect)),1);
end
if isempty(foci)
	nofocus=1;
	bestdep=initdep;
	foci=initfoci;
	bestcor=-1;bestrd=10;bestshift=0;
	return;
end

%% Determine the shifts for each focus ---------------------

if usecor
	mshift=mean(shiftCOR(foci));
	bestshift=shiftCOR(foci);
else
	mshift=mean(shiftRD(foci));
	bestshift=shiftRD(foci);	
end

bestcor=max(cors);
bestrd=min(rds(rds>0));

usecor=0;
isf=(0.5*length(foci))*fullfact(round(15/(0.5*length(foci)))*ones(size(foci))');
isf=isf+max(0,(mshift-max(max(isf))/2));
isf=isf-min(min(isf));
[a,b]=find(isf>maxshift);
isf(a,:)=[];

depi=DIST2W(:,foci);
for k=1:length(foci),
	depi(foci(k),k)=0;
end
for shifti=1:length(isf)
	if ~isempty(initdep)
		dep=min([initdep depi+ones(length(DIST2W),1)*isf(shifti,:)],[],2);
	else
		dep=min(depi+ones(length(DIST2W),1)*isf(shifti,:),[],2);
	end
    S=getSvplateau(T,dep,rep,GEOM.pS,[],scanmode);
    PSIA=lowpassma(AMA*S,lpass);
    usetime=min(qrsduration,min(minusetime,floor(qrsduration^2/max(dep))));        
    COR=corrcoef(PSIA(:,1:usetime),PSIREF(leads,1:usetime)); COR=COR(2,1);
	RD=norm(PSIREF(leads,1:usetime)-PSIA(:,1:usetime),'fro')/normphi;
	if COR > bestcor 
		bestcor=COR;
		if usecor
			bestshift=isf(shifti,:)';
% 			disp(['cor ' num2str([bestcor RD max(dep) foci' isf(shifti,:)],4)])
		end
	end
	if RD < bestrd
		bestrd=RD;
		if ~usecor
			bestshift=isf(shifti,:)';
			%             disp(['rd ' num2str([ bestrd COR max(dep) foci' max(dep) isf(shifti,:)],4)])
		end
	end
end

depi=DIST2W(:,foci);for k=1:length(foci),depi(foci(k),k)=0;end
if ~isempty(initdep)
	bestdep=min([initdep depi+ones(length(DIST2W),1)*bestshift'],[],2);
else
	bestdep=min(depi+ones(length(DIST2W),1)*bestshift',[],2);
end

% put new data in output variables

foci=[initfoci; foci];
bestshift=mean(bestshift);

%%
function [bestdep,bestshift]=optimizeShifts(GEOM,DIST2W,foci,shifts,leadset)
%% Determine the shifts for each focus ---------------------
scanmode=1;
global lpass

switch leadset
    case 'all';
        leads=1:size(GEOM.BSM,1);
    case'v12';
        leads=GEOM.v12;
    case'v12back'
        leads=GEOM.v12back;
end
onsetqrs  =GEOM.specs(2);
endqrs    =GEOM.specs(3);
qrsduration=endqrs-onsetqrs+1;
usetime=qrsduration;

% only p(1), setting the upslope of the TMP, is effective in this
% program (mode of gets=1)
t=0:qrsduration-1;
T=ones(length(GEOM.VER),1)*t;

orgdep=DIST2W(foci,:);
for k=1:length(foci),
	orgdep(k,:)=orgdep(k,:)+shifts(k);
end


PSIREF=baselinecor(lowpassma(GEOM.BSM(:,GEOM.specs(2):end),3));
PSIREF=PSIREF(:,1:qrsduration);
normphi=norm(PSIREF(leads,:),'fro');
AMA=zeromean(GEOM.AMA(leads,:));

%%
c=zeros(size(orgdep));
for i=1:length(foci)
	dep=orgdep;
	dep(i,:)=[];
	dep=min(dep);
    S=getSvplateau(T,dep,300*ones(size(dep)),GEOM.pS,[],scanmode);
    PSIA=lowpassma(AMA*S,lpass);
	COR=corrcoef(PSIA(:,1:usetime),PSIREF(leads,1:usetime)); 
	COR=COR(2,1);
	c(i)=COR;
end


%%
bestcor=-1;
bestrd=10;
usecor=1;
% shiftsi=4*ff2n(length(foci));
keep=ones(size(foci))'*0;
tkeep=ones(size(foci))'*0;
shifts=keep;
bestshift=shifts;
shiftsi=4*fullfact(round((qrsduration-70)/4)*ones(size(foci))');
for i=1:length(shiftsi)
	shifts=shiftsi(i,:);
	dep=min(DIST2W(:,foci)+ones(length(DIST2W),1)*shifts,[],2);
	S=getSvplateau(T,min(dep,[],2),300*ones(size(dep)),GEOM.pS,[],scanmode);
	PSIA=lowpassma(AMA*S,lpass);
	COR=corrcoef(PSIA(:,1:usetime),PSIREF(leads,1:usetime)); COR=COR(2,1);
	RD=norm(PSIREF(leads,1:usetime)-PSIA(:,1:usetime),'fro')/normphi;
	if COR > bestcor
		if COR > bestcor,bestcor=COR;end
		if usecor
			tkeep=keep|(shifts-i<=0);					
			bestshift=shifts;
			disp(['cor ' num2str([ bestcor RD max(dep) max(dep) shifts keep],4)])
		end
	end
	if RD < bestrd
		bestrd=RD;
		if ~usecor
			bestshift=shifts;
			tkeep=tkeep|shifts>i;
			disp(['rd ' num2str([ bestrd COR max(dep) max(dep) shifts keep],4)])        
		end
	end
	keep=tkeep;
end
foci(bestshift>qrsduration-80)=[];
bestshift(bestshift>qrsduration-80)=[];
bestdep=DIST2W(foci,:);
for k=1:length(foci),
	bestdep(k,:)=bestdep(k,:)+bestshift(k);
end
bestdep=min(bestdep)';

