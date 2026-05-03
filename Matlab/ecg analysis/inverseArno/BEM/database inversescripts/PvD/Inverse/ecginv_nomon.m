function [depfinal,repfinal,notchPotfinal,corfinal,rdfinal,iterfinal]=ecginv(GEOM,inittimes,initcase,runcase,usetimes,leads,mode)
% inverse determination of timing of Equivalent Double Layer source EDL; 
% optimization based on basic Marquardt type of solving thenon-linear 
% estimation problem alternate iterations between depolarization and 
% repolarization
% peter van Dam; 20080929

% if INV.mode==1, loopstart=1;nloop=1; repscore=0;end % dep only
% if INV.mode==2, loopstart=2;nloop=2; depscore=0;end % rep only
% if INV.mode==3, loopstart=3;nloop=3; end % notch
% if INV.mode==4, loopstart=1;nloop=2; end % dep rep 
% if INV.mode==5, loopstart=1;nloop=3; end % dep rep & notch
% if INV.mode==6, loopstart=6;nloop=6; end % AP amplitude only
global lpass	
global dodeprep

versie=1;	
casedir	='invresults\';

% PARAMETERS
lpass=5;				% moving avarage lowpass filtering
sigscal=GEOM.specs(1);				% for SIGPLOT
useSurfLapl=1;
maxiter=100;  
if strfind(GEOM.subject,'ppd1'),maxiter=10;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE INPUT DATA

% INV stucture
INV.VER=GEOM.VER;
INV.ITRI=GEOM.ITRI;

INV.mode=mode;
INV.BSM=GEOM.BSM(:,GEOM.specs(2):end);
INV.usetimes=min(usetimes,size(INV.BSM,2));
INV.BSM=baselinecor(INV.BSM,1,size(INV.BSM,2));
INV.ADJ=GEOM.ADJ;
if INV.mode==1,INV.BSM=baselinecor(INV.BSM,1,INV.usetimes);end
INV.PHI=zeromean(INV.BSM(leads,1:INV.usetimes));

INV.normphi=norm(INV.PHI(:,1:INV.usetimes),'fro');
INV.T=ones(size(GEOM.AMA,2),1)*(1:INV.usetimes);
% compute regularization operator, e.g.surface Laplacian
[INV.REGOP,INV.REGOPREP]=calcREGOP(GEOM,useSurfLapl);
%precompute large matrixes
INV.AMA=zeromean(GEOM.AMA(leads,:));
INV.ATA=INV.AMA'*INV.AMA;  
INV.REG=INV.REGOP'*INV.REGOP;
INV.REGREP=INV.REGOPREP'*INV.REGOPREP;

% criteria for ending iterations
INV.lambdamax=500;
INV.stopcrit=2e-4;
INV.pSn=GEOM.pS;

%% read initial estimates; 

% regularization parameters; must be used to tune the final result:
% small values may produce spatially 'wild' solutions: use trial and error
lambopt=0.1;  
% depolarization
OPT.DEP.mu=1.5e-4;   
if strcmp(GEOM.type,'_atria'),OPT.DEP.mu=1.e-6; end
OPT.DEP.lambopt=lambopt;  
OPT.DEP.tims=inittimes(:,1);
OPT.DEP.type='dep';

% repolarization
OPT.REP.tims=inittimes(:,2);
OPT.REP.lambopt=lambopt;  
OPT.REP.type='rep';
OPT.REP.mu=1.5e-4; 

% Jpoint
OPT.NOT.pot=zeros(size(GEOM.VER,1),1);
qrsduration=GEOM.specs(3)-GEOM.specs(2)+1;
if INV.mode==3 || INV.mode==5 || INV.mode==6
	OPT.NOT.pot=0.01*ones(size(GEOM.VER,1),1);
end
if dodeprep==0
	OPT.NOT.pot=1-sqrAmpl(GEOM,GEOM.specs(3)+10:GEOM.specs(3)+40,leads);
else
	OPT.AMP.pot=zeros(size(GEOM.VER,1),1);
end
OPT.NOT.usetimes=1:1:size(INV.PHI,2)-100;
OPT.NOT.lambopt=lambopt;
OPT.NOT.mu=1e-6;
OPT.NOT.type='notch';

% AP amplitude
if dodeprep==0
OPT.AMP.pot=ones(size(GEOM.VER,1),1);
else
OPT.AMP.pot=sqrAmpl(GEOM,GEOM.specs(3)+10:GEOM.specs(3)+40,leads);
end
if INV.mode==7,OPT.AMP.pot=OPT.AMP.pot-2;end
OPT.AMP.usetimes=1:1:size(INV.PHI,2)-100; %qrsduration+10:qrsduration+100;
OPT.AMP.lambopt=lambopt;
OPT.AMP.mu=1e-6;
OPT.AMP.type='amplitude';

%%
OPTSTART=OPT.DEP;
TST=gettres_v(INV,OPTSTART,OPT.REP,OPT.NOT,OPT.AMP);
%%    
% PRELUDE TO ITERATIVE OMPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute reg and res and tres for initial estimate

txt=sprintf('%s\n','iter   min   range  max    minAPD rangeAPD maxAPD reg    rd   ');
disp(txt);
temp=[0 min(OPT.DEP.tims) diff(range((OPT.DEP.tims))) max(OPT.DEP.tims) ...
		min(OPT.REP.tims-OPT.DEP.tims) diff(range(OPT.REP.tims-OPT.DEP.tims)) ...
		max(OPT.REP.tims-OPT.DEP.tims) TST.reg TST.rd ];
txt=sprintf('%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f Start time: %s',...
				temp,datestr(now,'HH,MM.SS'));
disp(txt);
 

%% 
% iterative approach to solving non-linear parameter estimation
% alternate between dep and rep

iter=0;score=1; % start outer loop
startTime=clock;
depscore=1;repscore=1;
if INV.mode==1, loopstart=1;nloop=1; repscore=0;end	% dep only
if INV.mode==2, loopstart=2;nloop=2; depscore=0;end % rep only
if INV.mode==3, loopstart=3;nloop=3; end			% notch only
if INV.mode==4, loopstart=1;nloop=2; end			% dep & rep 
if INV.mode==5, loopstart=1;nloop=2; end			% dep & notch & rep 
if INV.mode==6, loopstart=1;nloop=1; repscore=0;end % dep & notch
if INV.mode==7, loopstart=7;nloop=7; end			% AP amplitude only
totscore=0;
while iter<maxiter && (depscore || repscore) 
	iter=iter+1;
	for loop=loopstart:nloop
		
	switch loop
	case 1 % DEP
		[depscore,OPT.DEP,TST]=optimizeDepRep(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);
		if INV.mode==5 || INV.mode==6
			if dodeprep
			[tmpscore,OPT.AMP,TST]=optimizeAmpl(INV,OPT.AMP,OPT.DEP,OPT.REP,OPT.NOT);	
			else
			[tmpscore,OPT.NOT,TST]=optimizeNotch(INV,OPT.NOT,OPT.DEP,OPT.REP,OPT.AMP);				
			end
		end
		opt=OPT.DEP;		
		totscore=totscore+depscore;
	case 2 %REP
		if INV.mode==8
			[tmpscore,OPT.NOT]=optimizeNotch(INV,OPT.NOT,OPT.DEP,OPT.REP,OPT.AMP);				
		end		
		[repscore,OPT.REP,TST]=optimizeDepRep(INV,OPT.REP,OPT.DEP,OPT.NOT,OPT.AMP);
		opt=OPT.REP;
		totscore=totscore+repscore;				
	case 3
		[tmpscore,OPT.NOT,TST]=optimizeNotch(INV,OPT.NOT,OPT.DEP,OPT.REP,OPT.AMP);
		opt=OPT.DEP;															
	case 7
		[depscore,OPT.DEP,TST]=optimizeDepRep(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);
		[tmpscore,OPT.AMP,TST]=optimizeAmpl(INV,OPT.AMP,OPT.DEP,OPT.REP,OPT.NOT);	
		[repscore,OPT.REP,TST]=optimizeDepRep(INV,OPT.REP,OPT.DEP,OPT.NOT,OPT.AMP);				
		opt=OPT.DEP;		
		totscore=totscore+depscore;
	end
	temp=[iter min(opt.tims) diff(range(opt.tims)) max(opt.tims) ...
		  min(OPT.REP.tims-OPT.DEP.tims) diff(range(OPT.REP.tims-OPT.DEP.tims))...
		  max(OPT.REP.tims-OPT.DEP.tims) ...
		  TST.reg TST.rd ];
	txt=sprintf('%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f run %s time %s',...
				temp,opt.type,datestr(datenum(clock)-datenum(startTime),'HH,MM.SS'));
	disp(txt);
	end
end 
%%%% end outer loop; (loop1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TST=gettres_v(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);
depfinal=OPT.DEP.tims;
repfinal=OPT.REP.tims;
if dodeprep
	notchPotfinal=OPT.AMP.pot;
else
	notchPotfinal=OPT.NOT.pot;
end
rdfinal=TST.rd;
COR=corrcoef(TST.PHIA(:,1:INV.usetimes),INV.PHI(:,1:INV.usetimes));
corfinal=COR(2,1);
iterfinal=iter;
saveasci(outfil,[depfinal repfinal]);
saveasci(outfilP,INV.pSn);
fclose(fh);
% end of the inverse proceedure; view/and store results


%% =======================================================================
function [score,amplopt,TST]=optimizeAmpl(INV,amplopt,depopt,repopt,notchopt)

% prepare compution Marquardt step aimed at improving tims; 
% compute Sprime and Sprime*Sprime' based on previous iteration
TST=gettres_pot(INV,amplopt,depopt,repopt,notchopt);
starttres=TST.tres;
score=0;
lamb=amplopt.lambopt; 

Splus=getS(INV.T,depopt.tims,repopt.tims,INV.pSn,notchopt.pot,amplopt.pot+.001,4);
Smin =getS(INV.T,depopt.tims,repopt.tims,INV.pSn,notchopt.pot,amplopt.pot-.001,4);
Sprime=(Splus(:,amplopt.usetimes)-Smin(:,amplopt.usetimes))/2.;
SST=Sprime*Sprime';
% compute GTGM, with mu setting the weight of the regularization
% (smoothing based on the surface Laplacian)
GTGM=INV.ATA.*SST+amplopt.mu^2*INV.REG;
% compute gtres
gtres=sum(INV.AMA.*(TST.RES*Sprime'))';
righths=gtres-amplopt.mu^2*INV.REG*amplopt.pot;

% compute next estimate; test residual;
% increase lambda it until starttres < TST.tres
while TST.tres >= starttres, % inner loop
	lamb=lamb*2;
	if lamb > INV.lambdamax,  break, end
	GTGML=GTGM+lamb^2*eye(size(INV.AMA,2));
	deltau=inv(GTGML)*righths;
% 	weight function to limit plateau phase apmlitue to at least 85%
	newpot=amplopt.pot+deltau;
	newpot(newpot<0)=0;
	newpot(newpot>1)=1;
	testopt=amplopt;test.pot=newpot;
	TST=gettres_pot(INV,testopt,depopt,repopt,notchopt);
	if TST.tres<starttres
		if (starttres-TST.tres)/starttres >= INV.stopcrit
			ampl.pot=newpot; 
			score=1;
			disp('optimize amplitude')			
		end
		break;
	end
end  % end of inner loop
amplopt.lambopt=lamb/4;
TST=gettres_v(INV,depopt,repopt,notchopt,amplopt);
%% =======================================================================
function [score,notchopt,TST]=optimizeNotch(INV,notchopt,depopt,repopt,amplopt)

% prepare compution Marquardt step aimed at improving tims; 
% compute Sprime and Sprime*Sprime' based on previous iteration
TST=gettres_pot(INV,notchopt,depopt,repopt,amplopt);
starttres=TST.tres;
score=0;
lamb=notchopt.lambopt; 

Splus=getS(INV.T,depopt.tims,repopt.tims,INV.pSn,notchopt.pot+0.01,amplopt.pot,4);
Smin =getS(INV.T,depopt.tims,repopt.tims,INV.pSn,notchopt.pot-0.01,amplopt.pot,4);
Sprime=(Splus(:,notchopt.usetimes)-Smin(:,notchopt.usetimes))/2.;
SST=Sprime*Sprime';
% compute GTGM, with mu setting the weight of the regularization
% (smoothing based on the surface Laplacian)
GTGM=INV.ATA.*SST+notchopt.mu^2*INV.REG;
% compute gtres
gtres=sum(INV.AMA.*(TST.RES*Sprime'))';
righths=gtres-notchopt.mu^2*INV.REG*notchopt.pot;

% compute next estimate; test residual;
% increase lambda it until starttres < TST.tres

while TST.tres >= starttres, % inner loop
	lamb=lamb*2;
	if lamb > INV.lambdamax,  break, end
	GTGML=GTGM+lamb^2*eye(size(INV.AMA,2));
	delnotch=inv(GTGML)*righths;
% 	weight function to limit plateau phase apmlitue to at least 85%
	newpot=notchopt.pot+delnotch; 
% 10 is appoximates the (liniear?) relation between notchpotential 
% 	and the parameter regulating the notch amplitude
	
	newpot(newpot<0)=0;
	newpot(newpot>1)=1;
	testopt=notchopt;testopt.pot=newpot;
	TST=gettres_pot(INV,testopt,depopt,repopt,amplopt);
	if TST.tres<starttres
		if (starttres-TST.tres)/starttres >= INV.stopcrit
			notchopt.pot=newpot; 
			score=1;
			disp('optimize notch')
		end
		break;
	end
end  % end of inner loop
notchopt.lambopt=lamb/4;
TST=gettres_v(INV,depopt,repopt,notchopt,amplopt);

%% =======================================================================
function [score,opt,TST]=optimizeDepRep(INV,opt,keepopt,notchopt,amplopt)

% prepare compution Marquardt step aimed at improving tims; 
% compute Sprime and Sprime*Sprime' based on previous iteration

TST=gettres_v(INV,opt,keepopt,notchopt,amplopt);
starttres=TST.tres;
score=0;
lamb=opt.lambopt; 

if strcmp(opt.type,'dep')
Splus=getS(INV.T,opt.tims+1,keepopt.tims,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
Smin =getS(INV.T,opt.tims-1,keepopt.tims,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
startopt=opt;	
elseif strcmp(opt.type,'rep')
Splus=getS(INV.T,keepopt.tims,opt.tims+1,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
Smin =getS(INV.T,keepopt.tims,opt.tims-1,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
startopt=opt;
else %if strcmp(opt.type,'apd')
Splus=getS(INV.T,keepopt.tims,opt.tims+1,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
Smin =getS(INV.T,keepopt.tims,opt.tims-1,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
startopt=opt;
startopt.tims=startopt.tims-keepopt.tims;
end
   
Sprime=(Splus(:,1:INV.usetimes)-Smin(:,1:INV.usetimes))/2.;
SST=Sprime*Sprime';
% compute GTGM, with mu setting the weight of the regularization
% (smoothing based on the surface Laplacian)
if strcmp(startopt.type,'rep')
	GTGM=INV.ATA.*SST+opt.mu^2*INV.REGREP;
	% compute gtres
	gtres=sum(INV.AMA.*(TST.RES*Sprime'))';
	righths=gtres-opt.mu^2*INV.REGREP*startopt.tims;	
else
	GTGM=INV.ATA.*SST+opt.mu^2*INV.REG;
	% compute gtres
	gtres=sum(INV.AMA.*(TST.RES*Sprime'))';
	righths=gtres-opt.mu^2*INV.REG*startopt.tims;
end
testopt=startopt;
% compute next estimate; test residual;
% increase lambda it until starttres < TST.tres

while TST.tres >= starttres, % inner loop
	lamb=lamb*2;
	if lamb > INV.lambdamax,  break, end
	GTGML=GTGM+lamb^2*eye(size(INV.AMA,2));
	deltau=pinv(GTGML)*righths;
	testopt.tims=startopt.tims+deltau; 
	if strcmp(startopt.type,'apd')
		testopt.tims=testopt.tims+keepopt.tims;
		TST=gettres_v(INV,testopt,keepopt,notchopt,amplopt);
	else
		TST=gettres_v(INV,testopt,keepopt,notchopt,amplopt);		
	end
	if TST.tres<starttres
		if (starttres-TST.tres)/starttres >= INV.stopcrit
			opt.tims=testopt.tims; 
			score=1;
		end
		break;
	end
end  % end of inner loop
opt.lambopt=lamb/4;

%% =======================================================================
function TST=gettres_pot(INV,opt,depopt,repopt,otheropt)

global lpass
if strcmp(opt.type,'notch')
	S=getS(INV.T,depopt.tims,repopt.tims,INV.pSn,opt.pot,otheropt.pot,INV.mode);
else %type == ' ampl',
	S=getS(INV.T,depopt.tims,repopt.tims,INV.pSn,otheropt.pot,opt.pot,INV.mode);
end
TST.PHIA=baselinecor(lowpassma(INV.AMA*S(:,opt.usetimes),lpass));
TST.RES=INV.PHI(:,opt.usetimes)-TST.PHIA;
TST.rd=norm(TST.RES,'fro')/norm(INV.PHI(:,opt.usetimes),'fro'); %NOTE: unfiltered rd
TST.reg=norm(INV.REGOP*opt.pot)/1000;
TST.tres=sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
TST.S=S;

%% =======================================================================
function TST=gettres_v(INV,opt,keepopt,notchopt,amplopt)

global lpass
if strcmp(opt.type,'dep')
	S=getS(INV.T,opt.tims,keepopt.tims,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
else %pol=='rep',
	S=getS(INV.T,keepopt.tims,opt.tims,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
end
TST.PHIA=baselinecor(lowpassma(INV.AMA*S(:,1:INV.usetimes),lpass));
TST.RES=INV.PHI(:,1:INV.usetimes)-TST.PHIA;
TST.rd=norm(TST.RES,'fro')/INV.normphi; %NOTE: unfiltered rd
if strcmp(opt.type,'rep')
	TST.reg=norm(INV.REGOPREP*opt.tims)/1000;
else
	TST.reg=norm(INV.REGOP*opt.tims)/1000;
end
TST.tres=sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
TST.S=S;

%% =======================================================================
function [REGOP,REGOPREP]=calcREGOP(GEOM,useSurfLapl)

if useSurfLapl==1
	if max(max(abs(GEOM.VER)))>1
		REGOP=surflapl(GEOM.VER/1000,GEOM.ITRI,1);
	else
		REGOP=surflapl(GEOM.VER,GEOM.ITRI,1);
	end
	REGOPREP=REGOP;
else %if useSurfLapl==2
	disp('inverse 1.distance^2 regularization');
% 	REGOP=GEOM.DIST2W;
%     REGOP=GEOM.ANISDIST;
%     ADJ=GEOM.ADJ; ADJ(ADJ>50)=0;
%     REGOP=calcAnisADJ(GEOM.VER,GEOM.ITRI,ADJ,GEOM.RVER,GEOM.LVER,2.5);
    REGOP=calcAnisADJ(GEOM,2.5);
	if max(max(abs(GEOM.VER)))>1
		REGOP=REGOP/1000;
	end
	REGOP=REGOP.^(-2);
	REGOP(REGOP==Inf)=0;
	scal=sum(REGOP,2);
	scal=1./scal;
	REGOP=diag(scal)*REGOP;
    % 59.7 scales the result to be of the same order of SL
	REGOP=70*(REGOP-eye(size(REGOP)));% 
end

%% =======================================================================
   
function S=getS(T,dep,rep,p,notchpot,scaleAmpl,mode)
global dodeprep
if dodeprep
S=getSmode(T,dep,rep,p,scaleAmpl,mode);
else
S=getSmode(T,dep,rep,p,notchpot,mode);
end