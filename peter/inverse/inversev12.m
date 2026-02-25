
function meas=inversev12(varargin)%GEOM,inittimes,initcase,runcase,usetimes,leads,mode,mudep,murep)
% inverse determination of timing of Equivalent Double Layer source EDL; 
% optimization based on basic Marquardt type of solving thenon-linear 
% estimation problem alternate iterations between depolarization and 
% repolarization
% Peter van Dam; 2010 november. 
% All rights reserved Peacs, Arnhem  the Netherlands

% if INV.mode==1, loopstart=1;nloop=1; repscore=0;end % dep only
% if INV.mode==2, loopstart=2;nloop=2; depscore=0;end % rep only
% if INV.mode==3, loopstart=3;nloop=3; end % notch
% if INV.mode==4, loopstart=1;nloop=2; end % dep rep 
% if INV.mode==5, loopstart=1;nloop=3; end % dep rep & notch
% if INV.mode==6, loopstart=6;nloop=6; end % AP amplitude only
global lpass	
global doAmplitude
doAmplitude = 0;
% inverse procedure parameters
mudep = 1.5e-4;
murep = 1.5e-4;
muampl = 1.5e-4;

repOpt='apd';
MINRD=0.15;
maxiter=25; 
lambopt=0.1;  
useAmpl = 0;
useNotch =0;
casedir	='invresults\';
subname='';
leads = 1:64;
mode = 4;
doWeight = 0;
%%
if length(varargin) < 3
	error('This routine needs at least three parameters');
else
    GEOM = varargin{1};
    initdep = varargin{2};
    initrep = varargin{3};
    usetimes = size(GEOM.BSM,2)-GEOM.specs(2)+1;
	pp=4;
    while pp<=nargin
        if ischar(varargin{pp})
            key=lower(varargin{pp});
            switch key
                case 'mudep'
                    mudep = varargin{pp+1};pp=pp+2;
                case 'murep'
                    murep = varargin{pp+1};pp=pp+2;
                case 'repopt'
                    repOpt = varargin{pp+1};pp=pp+2;                    
                case 'muampl'
                    muampl = varargin{pp+1};pp=pp+2;                    
                case 'estimateampl'
                    useAmpl = varargin{pp+1};pp=pp+2;
                case 'estimatenotch'
                    useNotch = varargin{pp+1};pp=pp+2;
                case 'maxiter'
                    maxiter = varargin{pp+1};pp=pp+2;
                case 'minrd'
                    MINRD = varargin{pp+1};pp=pp+2;
                case 'casedir'
                    casedir = varargin{pp+1};pp=pp+2;			
                case 'logname'
                    subname = varargin{pp+1};pp=pp+2;
                case 'usetimes'
                    usetimes = varargin{pp+1};pp=pp+2;                    
                case 'leads'
                    leads = varargin{pp+1};pp=pp+2;                    
                case 'weighed'
                    doWeight = varargin{pp+1};pp=pp+2;  
                case 'mode'
                    mode = varargin{pp+1};pp=pp+2;  
                otherwise
                    error('unknown parameter');
            end
        end
    end
end

%%
if ~strcmp(repOpt,'rep')
    repOpt = 'apd';
end

if ~exist(casedir,'dir')
    mkdir(casedir);
end
runcase = [casedir GEOM.subject GEOM.beat subname];
versie=1;	


% PARAMETERS
lpass=5;				% moving avarage lowpass filtering
useSurfLapl=1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE INPUT DATA

% INV stucture
INV.VER=GEOM.VER;
INV.ITRI=GEOM.ITRI;

INV.mode=mode;
INV.foci = calcLocmin(GEOM.DIST,initdep)';


if strcmp(GEOM.type,'_atria') % no baseline coroection for the atria
    INV.BSM=(GEOM.BSM(:,GEOM.specs(2):end));
    if INV.mode == 1
         INV.usetimes=min(GEOM.specs(3)-GEOM.specs(2)+1,size(INV.BSM,2));
    else
        INV.usetimes=min(usetimes,size(INV.BSM,2));
    end
    INV.PHIREF=(INV.BSM(leads,1:INV.usetimes));
else
    INV.BSM=baselinecor((GEOM.BSM(:,GEOM.specs(2):end)));
    if INV.mode == 1
        INV.usetimes=min(GEOM.specs(3)-GEOM.specs(2)+1,size(INV.BSM,2));
    else
        INV.usetimes=min(usetimes,size(INV.BSM,2));
    end
    INV.PHIREF=baselinecor((INV.BSM(leads,1:INV.usetimes)));
end

INV.ADJ=GEOM.ADJ;
INV.ADJsurf = GEOM.ADJsurf;
INV.ADJ2W=GEOM.ADJ2W;
INV.DIST=GEOM.DIST;

INV.normphi=norm(INV.PHIREF,'fro');
INV.T=ones(size(GEOM.AMA,2),1)*(0:INV.usetimes-1);

% compute regularization operator, e.g.surface Laplacian
[INV.REGOP,INV.REGOPREP]=calcREGOP(GEOM,useSurfLapl);
%precompute large matrixes
INV.AMA=(GEOM.AMA(leads,:));
INV.ATA=INV.AMA'*INV.AMA;  

if doWeight
   weight = rms(INV.PHIREF')';
    weight = weight/ max(weight);
    weight( weight < max(weight)/1.5) = max(weight)/1.5;
    INV.PHIREF = INV.PHIREF .* (weight * ones(1,size(INV.PHIREF,2)));
    INV.AMA = INV.AMA .*(weight * ones(1,size(INV.AMA,2)));
end
INV.REG = INV.REGOP'*INV.REGOP;
INV.REGREP = INV.REGOPREP'*INV.REGOPREP;


INV.pSn=GEOM.pS;
% criteria for ending iterations
INV.lambdamax=500;
INV.stopcrit=2e-4;
%% read initial estimates; 

% regularization parameters; must be used to tune the final result:
% small values may produce spatially 'wild' solutions: use trial and error

% depolarization
OPT.DEP.mu=mudep; %1.4e-4  
OPT.DEP.lambopt=lambopt;  
OPT.DEP.tims=initdep;
OPT.DEP.type='dep';

% repolarization
OPT.REP.tims=initrep;
OPT.REP.lambopt=lambopt;  
OPT.REP.type=repOpt;
OPT.REP.mu=murep; 

% Jpoint
OPT.NOT.pot=zeros(size(GEOM.VER,1),1);
if useNotch% future when also an notch is included in the source model
    OPT.NOT.pot=1-sqrAmpl(GEOM,GEOM.specs(3)+10:GEOM.specs(3)+40,leads);
    OPT.NOT.lambopt=lambopt;
    OPT.NOT.mu=1e-6;
    OPT.NOT.type='notch';
end
qrsduration=GEOM.specs(3)-GEOM.specs(2)+1;

% AP amplitude
OPT.AMP.pot = ones(size(GEOM.VER,1),1);
if useAmpl
    OPT.AMP.pot = sqrAmpl(GEOM,GEOM.specs(3)+10:GEOM.specs(3)+40,leads);
    OPT.AMP.usetimes=qrsduration+1:qrsduration+60;
    OPT.AMP.lambopt=lambopt;
    OPT.AMP.mu=muampl;
    OPT.AMP.type='amplitude';
end
%%
OPTSTART=OPT.DEP;
TST=gettres_v(INV,OPTSTART,OPT.REP,OPT.NOT,OPT.AMP);
%% =======================================================================
% DOCUMENT CASE  
saveasci([runcase 'srcinit'],[initdep initrep]);
logfile=[runcase '.log'];
outfil= [runcase '.src'];
outfilP= [runcase '.pS'];
fh=fopen(logfile,'w');
fprintf(fh,'%s\n',['file: ' logfile]);
fprintf(fh,'%s\n',['date: ' datestr(clock)]);
fprintf(fh,'%s\n',['prog: invgeom (' num2str(versie) ')' ]);
fprintf(fh,'%s\n',['initvals:   ' [runcase 'srcinit']]);
fprintf(fh,'%s\n',['subject:   ' GEOM.subject]);
fprintf(fh,'%s\n',['used anisotropyRatio:    ' GEOM.anisotropyRatio]); 

if useSurfLapl==1
	fprintf(fh,'%s\n','SL regularization');
else 
	fprintf(fh,'%s\n','inverse 1.distance^2 regularization ');
end
fprintf(fh,'mudep: %0.4g    murep: %0.4g \n',OPT.DEP.mu,OPT.REP.mu);
fprintf(fh,'%s\n','ifix(1) tfix(1) dampd   dampr usetimes   tblc   mode  maxiter');

RESNOW=[];
%%    
% PRELUDE TO ITERATIVE OMPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute reg and res and tres for initial estimate

fprintf(fh,'%s\n','iter   min    mean     max     STD    reg     rd     tresd');
temp=[0 min(OPT.DEP.tims) diff(range(OPT.DEP.tims)) max(OPT.DEP.tims) std(OPT.DEP.tims) TST.reg TST.rd TST.tres];
sprintf(   '%3d %7.1f %7.1f %7.1f  %6.1f %6.1f %7.4f %7.4f',temp(1:8));
fprintf(fh,'%3d %7.1f %7.1f %7.1f  %6.1f %6.1f %7.4f %7.4f\n',temp(1:8));

fprintf(fh,'%s\n','iter   min   range  max    std  minAPD rangeAPD maxAPD stdAPD  reg    rd     tresd');
txt=sprintf('%s\n','iter   min   range  max    STD  minAPD rangeAPD maxAPD stdAPD  reg    rd     tresd');
disp(txt);
temp=[0 min(OPT.DEP.tims) diff(range((OPT.DEP.tims))) max(OPT.DEP.tims) std(OPT.DEP.tims)...
				 min(OPT.REP.tims-OPT.DEP.tims) diff(range(OPT.REP.tims-OPT.DEP.tims)) ...
				 max(OPT.REP.tims-OPT.DEP.tims) std(OPT.REP.tims-OPT.DEP.tims)...
				 TST.reg TST.rd TST.tres];
txt=sprintf(   '%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f Start time: %s',temp,datestr(now,'HH,MM.SS'));
disp(txt);
fprintf(fh,'%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f',temp);
fprintf(fh,'time run: %s\n',datestr(now,'HH,MM.SS')) ;
RESNOW(1,:)=[temp norm(TST.RES)];
ik=2;
 

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
TST.rd=1;
while iter<maxiter && (depscore || repscore) && TST.rd>MINRD 
	iter=iter+1;
	for loop=loopstart:nloop
		switch loop
			case 1 % DEP
				[depscore,OPT.DEP,TST]=optimizeDepRep(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);
				if INV.mode==5 || INV.mode==6
					if doAmplitude
					[tmpscore,OPT.AMP,TST]=optimizeAmpl(INV,OPT.AMP,OPT.DEP,OPT.REP,OPT.NOT);	
					else
					[tmpscore,OPT.NOT,TST]=optimizeNotch(INV,OPT.NOT,OPT.DEP,OPT.REP,OPT.AMP);				
					end
				end
				opt=OPT.DEP;		
                score = depscore;
				totscore=totscore+depscore;
			case 2 %REP
				if INV.mode==8
					[tmpscore,OPT.NOT,TST]=optimizeNotch(INV,OPT.NOT,OPT.DEP,OPT.REP,OPT.AMP);				
				end		
				[repscore,OPT.REP,TST]=optimizeDepRep(INV,OPT.REP,OPT.DEP,OPT.NOT,OPT.AMP);
				opt=OPT.REP;
                score = repscore;
				totscore=totscore+repscore;				
			case 3
				[tmpscore,OPT.NOT,TST]=optimizeNotch(INV,OPT.NOT,OPT.DEP,OPT.REP,OPT.AMP);
				opt=OPT.DEP;															
			case 7
				[depscore,OPT.DEP,TST]=optimizeDepRep(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);
				[tmpscore,OPT.AMP,TST]=optimizeAmpl(INV,OPT.AMP,OPT.DEP,OPT.REP,OPT.NOT);	
				[repscore,OPT.REP,TST]=optimizeDepRep(INV,OPT.REP,OPT.DEP,OPT.NOT,OPT.AMP);				
				opt=OPT.DEP;		
                score = depscore;
				totscore=totscore+depscore;
		end
		temp=[iter min(opt.tims) diff(range(opt.tims)) max(opt.tims) std(opt.tims)...
			  min(OPT.REP.tims-OPT.DEP.tims) diff(range(OPT.REP.tims-OPT.DEP.tims))...
			  max(OPT.REP.tims-OPT.DEP.tims) std(OPT.REP.tims-OPT.DEP.tims)...
			  TST.reg TST.rd TST.tres];
		if score==0
			k=rem(ik,2)+1;
			temp(9:10)=RESNOW(ik-k,9:10);
			temp(11)=RESNOW(ik-1,11);
		end
		txt=sprintf(   '%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f run %s time %s',...
					temp,opt.type,datestr(datenum(clock)-datenum(startTime),'HH,MM.SS'));
		disp(txt);
		fprintf(fh,'%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f',temp);
		fprintf(fh,'run %s time %s\n',opt.type,datestr(datenum(clock)-datenum(startTime),'HH,MM.SS'));
		RESNOW(ik,:)=[temp norm(TST.RES)];
		ik=ik+1;
	end
end 
%%%% end outer loop; (loop1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TST=gettres_v(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);
meas.depfinal=OPT.DEP.tims;
meas.repfinal=OPT.REP.tims;
if useAmpl
	meas.amplfinal=OPT.AMP.pot;
end
if useNotch
	meas.notchPotfinal=OPT.NOT.pot;
end
meas.rdfinal=TST.rd;
COR=corrcoef(TST.PHIA,INV.PHIREF);
meas.corfinal=COR(2,1);
meas.iterfinal=iter;
meas.log=RESNOW;

saveasci(outfil,[meas.depfinal meas.repfinal]);
saveasci(outfilP,INV.pSn);
fclose(fh);


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

Sprime=(Splus-Smin)/2.;
SST=Sprime*Sprime';
% compute GTGM, with mu setting the weight of the regularization
% (smoothing based on the surface Laplacian)
if strcmp(startopt.type,'dep')
	GTGM=INV.ATA.*SST+opt.mu^2*INV.REG;
	% compute gtres
	gtres=sum(INV.AMA.*(TST.RES*Sprime'))';
	righths=gtres-opt.mu^2*INV.REG*startopt.tims;
else
    GTGM=INV.ATA.*SST+opt.mu^2*INV.REGREP;
	% compute gtres
	gtres=sum(INV.AMA.*(TST.RES*Sprime'))';
	righths=gtres-opt.mu^2*INV.REGREP*startopt.tims;	    
end
testopt=startopt;
% compute next estimate; test residual;
% increase lambda it until starttres < TST.tres

while TST.tres >= starttres, % inner loop
	lamb=lamb*2;
	if lamb > INV.lambdamax,  break, end
	GTGML=GTGM+lamb^2*eye(size(INV.AMA,2));
	deltau=pinv(GTGML)*righths;
	newtime=startopt.tims+deltau; 
	if strcmp(startopt.type,'dep')        
%         INV.foci = calcLocmin(INV.DIST,newtime,INV.foci);
        velo = calcActVelo(newtime,INV.ADJ2W,INV.DIST,INV.foci);
        foci =[];
        redo=false;
        for  k=1:length(velo)
            if velo(k) > 1 && any(INV.ADJsurf(k,INV.foci) > 0) 
                foci = [foci; k];
                redo=true;
            end
        end
        if redo
            INV.foci = unique([foci; INV.foci]);
            velo = calcActVelo(newtime,INV.ADJ2W,INV.DIST,INV.foci);
        end
        velo(INV.foci) = 0.5; % prevent adaptation
        deltau(velo > 0.6) = deltau(velo > 0.6 )./((velo(velo>0.6) - 0.55)*15);
        deltau(velo < 0.1 ) = 0;%-deltau(velo < 0.1)/2;
        deltau(velo > 1.0) = 0;%-deltau(velo > 1.5)/2;
        newtime = startopt.tims + deltau; 
		newtime(newtime<0)=0;        
	elseif strcmp(startopt.type,'rep')
		newtime(newtime<max(keepopt.tims)+50)=max(keepopt.tims)+50;
	elseif strcmp(startopt.type,'apd')
		newtime(newtime<50)=50;
	end	
	
	testopt.tims=newtime;
	if strcmp(startopt.type,'apd')
		testopt.tims=testopt.tims+keepopt.tims;
		TST=gettres_v(INV,testopt,keepopt,notchopt,amplopt);
	else
		TST=gettres_v(INV,testopt,keepopt,notchopt,amplopt);		
	end
	if TST.tres < starttres || sum(deltau) == 0
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
TST.PHIA=lowpassma(INV.AMA*S(:,opt.usetimes),lpass);
% TST.PHIA=baselinecor(lowpassma(INV.AMA*S(:,opt.usetimes),lpass));
TST.RES=INV.PHIREF(:,opt.usetimes)-TST.PHIA;
TST.rd=norm(TST.RES,'fro')/norm(INV.PHIREF(:,opt.usetimes),'fro'); %NOTE: unfiltered rd
TST.reg=norm(INV.REGOP*opt.pot)/1000;
TST.tres=sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
TST.S=S;

%% =======================================================================
function TST=gettres_v(INV,opt,keepopt,notchopt,amplopt)

global lpass
if strcmp(opt.type,'dep')
	TST.S=getS(INV.T,opt.tims,keepopt.tims,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
else %pol=='rep',
	TST.S=getS(INV.T,keepopt.tims,opt.tims,INV.pSn,notchopt.pot,amplopt.pot,INV.mode);
end
% TST.PHIA=baselinecor(lowpassma(INV.AMA*TST.S,lpass));
TST.PHIA=lowpassma(INV.AMA*TST.S,lpass);
TST.RES=INV.PHIREF-TST.PHIA(1:size(INV.PHIREF,1),1:size(INV.PHIREF,2));
TST.rd=norm(TST.RES,'fro')/INV.normphi; %NOTE: unfiltered rd
if strcmp(opt.type,'dep')
	TST.reg=norm(INV.REGOP*opt.tims)/1000;
else
	TST.reg=norm(INV.REGOPREP*opt.tims)/1000;
end
TST.tres=sqrt(TST.rd^2+(TST.reg*opt.mu)^2);

%% =======================================================================
function [REGOP,REGOPREP]=calcREGOP(GEOM,useSurfLapl)

% if isfield(GEOM,'LAPL')
%     REGOP2 = 2.0 * GEOM.LAPL;
%     REGOPREP2 = REGOP;
% end
if useSurfLapl==1
	if max(max(abs(GEOM.VER)))>1
		REGOP=surflapl(GEOM.VER/1000,GEOM.ITRI,1);
	else
		REGOP=surflapl(GEOM.VER,GEOM.ITRI,1);
    end
    if isfield(GEOM,'LAPL')
        REGOP = REGOP - 2 * GEOM.LAPL;
    end   
    LD =GEOM.ADJ2W;
    LD(LD > 0) = (LD( LD > 0 ).^-1.5);
    REGOPREP = LD - diag(sum(LD));


else %if useSurfLapl==2
	disp('inverse 1.distance^2 regularization');
%     LD =GEOM.ADJ2W;
    LD = GEOM.DIST;
    LD(LD > 50) = 0;
    LD(LD > 0) = (LD( LD > 0 ).^-1.5);
    REGOPREP = LD - diag(sum(LD));
    LD = GEOM.DIST;
    LD(LD > 30) = 0;
    LD(LD > 0) = (LD( LD > 0 ).^-1.5);   
    REGOP = LD;

end

%% =======================================================================
   
function S=getS(T,dep,rep,p,notchpot,scaleAmpl,mode)

if mode == 4
    if size(T,2) < max(rep) + 50
        Tt=ones(length(dep),1)*([0:max(rep) + 50]);     
        S=getSmode(Tt,dep,rep,p,scaleAmpl,4);
        S=S(:,1:size(T,2));
    else
        S=getSmode(T,dep,rep,p,scaleAmpl,4);
    end
else
    S=getSmode(T,dep,rep,p,notchpot,mode);
end