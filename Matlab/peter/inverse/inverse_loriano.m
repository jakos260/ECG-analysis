function meas=inverse_loriano(varargin)%GEOM,inittimes,initcase,runcase,usetimes,leads,mode,mudep,murep)
% inverse determination of timing of Equivalent Double Layer source EDL; 
% optimization based on basic Marquardt type of solving thenon-linear 
% estimation problem alternate iterations between depolarization and 
% repolarization
% Peter van Dam; 2010 november. 
% All rights reserved Peacs, Arnhem  the Netherlands

% if INV.mode==1, loopstart=1;nloop=1; repscore=0;end % dep only
% if INV.mode==2, loopstart=2;nloop=2; depscore=0;end % rep only
% if INV.mode==4, loopstart=1;nloop=2; end % dep rep 
global lpass	
% inverse procedure parameters
mudep = 1e-5;
murep = 1e-5;

repOpt='apd';
MINRD=0.15;
maxiter=25; 
lambopt=0.1;  
useAmpl = 0;
useNotch =0;
casedir	='invresults\';
subname='';
leads = 1:65;
mode = 4; % rep only
doWeight = 0;
reg=[];
%%
if length(varargin) < 3
	error('This routine needs at least three parameters');
else
    GEOM.VER  = varargin{1}; % heart vertices
    GEOM.ITRI = varargin{2}; % heart triangles
    GEOM.BSM  = varargin{3}; % BSM of the relevant part of the ECG (QRST)
    
    initdep = varargin{4}; % the original dep times
    initrep = varargin{5}; % the rep times
    GEOM.subject = varargin{6};  
    GEOM.AMA = varargin{7};
    GEOM.qrsduration = varargin{8};
    GEOM.pS = varargin{9}; % your tmp repolarization paramters
    
    usetimes = size(GEOM.BSM,2);
    
    INV.useWeighedRd = 1;
	pp=10;
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
                case 'reg'
                    reg = varargin{pp+1};pp=pp+2;  
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
subname = [ '_' num2str(length(leads))];
runcase = fullfile(casedir,[ GEOM.subject subname]); 
versie=1;	


% PARAMETERS
lpass=5;				% moving avarage lowpass filtering
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE INPUT DATA

% INV stucture
INV.VER=GEOM.VER;
INV.ITRI=GEOM.ITRI;

INV.mode=mode;
INV.BSM = GEOM.BSM;
INV.usetimes = usetimes;
INV.PHIREF = INV.BSM(:,1:INV.usetimes);

INV.normphi = norm(INV.PHIREF,'fro');
INV.T = ones(size(GEOM.AMA,2),1)*(0:INV.usetimes-1);

% compute regularization operator, e.g.surface Laplacian
if isempty(reg)
    [INV.REGOP,INV.REGOPREP]=calcREGOP(GEOM);
else
    INV.REGOP = reg;
    INV.REGOPREP = INV.REGOP;
end

INV.ROTRO = INV.REGOP'*INV.REGOP;
INV.ROTROREP = INV.REGOPREP'*INV.REGOPREP;

%precompute large matrixes
INV.AMA=zeromean(GEOM.AMA);
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
OPT.REP.type = 'apd';



%%
OPTSTART=OPT.DEP;
TST=gettres_v(INV,OPTSTART,OPT.REP);
%% =======================================================================
% DOCUMENT CASE  
saveasci([runcase '.srcinit'],[initdep initrep]);
logfile=[runcase '.log'];
outfil= [runcase '.src'];
fh=fopen(logfile,'w');
fprintf(fh,'%s\n',['file: ' logfile]);
fprintf(fh,'%s\n',['date: ' datestr(clock)]);
fprintf(fh,'%s\n',['prog: invgeom (' num2str(versie) ')' ]);
fprintf(fh,'%s\n',['initvals:   ' [runcase 'srcinit']]);
fprintf(fh,'%s\n',['subject:   ' GEOM.subject]);


fprintf(fh,'%s\n','SL regularization');
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
if INV.mode==4, loopstart=1;nloop=2; end			% dep & rep 
totscore=0;
TST.rd=1;
while iter<maxiter && (depscore || repscore) && TST.rd>MINRD 
	iter=iter+1;
	for loop=loopstart:nloop
		switch loop
			case 1 % DEP
				[depscore,OPT.DEP,TST]=optimizeDepRep(INV,OPT.DEP,OPT.REP);
				opt=OPT.DEP;		
                score = depscore;
				totscore=totscore+depscore;
			case 2 %REP
				[repscore,OPT.REP,TST]=optimizeDepRep(INV,OPT.REP,OPT.DEP);
				opt=OPT.REP;
                score = repscore;
				totscore=totscore+repscore;				
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
TST=gettres_v(INV,OPT.DEP,OPT.REP);
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

fclose(fh);

%% =======================================================================
function [score,opt,TST]=optimizeDepRep(INV,opt,keepopt)

% prepare compution Marquardt step aimed at improving tims; 
% compute Sprime and Sprime*Sprime' based on previous iteration

TST=gettres_v(INV,opt,keepopt);
starttres=TST.tres;
score=0;
lamb=opt.lambopt; 

if strcmp(opt.type,'dep')
	Splus=getS(INV.T,opt.tims+1,keepopt.tims,INV.pSn,INV.mode);
	Smin =getS(INV.T,opt.tims-1,keepopt.tims,INV.pSn,INV.mode);
	startopt=opt;	
elseif strcmp(opt.type,'rep')
	Splus=getS(INV.T,keepopt.tims,opt.tims+1,INV.pSn,INV.mode);
	Smin =getS(INV.T,keepopt.tims,opt.tims-1,INV.pSn,INV.mode);
	startopt=opt;
else %if strcmp(opt.type,'apd')
	Splus=getS(INV.T,keepopt.tims,opt.tims+1,INV.pSn,INV.mode);
	Smin =getS(INV.T,keepopt.tims,opt.tims-1,INV.pSn,INV.mode);
	startopt=opt;
	startopt.tims=startopt.tims-keepopt.tims;
end

Sprime=(Splus-Smin)/2.;
SST=Sprime*Sprime';
% compute GTGM, with mu setting the weight of the regularization
% (smoothing based on the surface Laplacian)
if strcmp(startopt.type,'dep')
    muREG = bsxfun(@times,opt.mu^2,INV.ROTRO);
    M = INV.ATA.*SST;
else
    muREG = bsxfun(@times,opt.mu^2,INV.ROTROREP);
    M = INV.ATA.*SST;
end
GTGM = M'*M + muREG;
% GTGM = M + muREG; old
gtres=sum(INV.AMA'.*(Sprime * TST.RES'),2);
righths= M' * gtres - muREG * startopt.tims;

testopt=startopt;
% compute next estimate; test residual;
% increase lambda it until starttres < TST.tres

while TST.tres >= starttres, % inner loop
	lamb=lamb*2;
	if lamb > INV.lambdamax,  break, end
	GTGML = GTGM + lamb^2*eye(size(INV.AMA,2));

    deltau=GTGML \ righths;  % INVERSE; total computing time  0.87 times preceding

    newtime=startopt.tims+deltau;

	newtime=startopt.tims+deltau; 
	if strcmp(startopt.type,'dep')
		newtime(newtime<-5)=-5;
	elseif strcmp(startopt.type,'rep')
		newtime(newtime<max(keepopt.tims)+50)=max(keepopt.tims)+50;
	elseif strcmp(startopt.type,'apd')
		newtime(newtime<50)=50;
	end	
	
	testopt.tims=newtime;
	if strcmp(startopt.type,'apd')
		testopt.tims=testopt.tims+keepopt.tims;
		TST=gettres_v(INV,testopt,keepopt);
	else
		TST=gettres_v(INV,testopt,keepopt);		
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
function TST=gettres_v(INV,opt,keepopt)

global lpass
if strcmp(opt.type,'dep')
	TST.S=getS(INV.T,opt.tims,keepopt.tims,INV.pSn,INV.mode);
else %pol=='rep',
	TST.S=getS(INV.T,keepopt.tims,opt.tims,INV.pSn,INV.mode);
end
% TST.PHIA=baselinecor(lowpassma(INV.AMA*TST.S,lpass));
TST.PHIA=lowpassma(INV.AMA*TST.S,lpass);
TST.RES=INV.PHIREF - TST.PHIA(1:size(INV.PHIREF,1),1:size(INV.PHIREF,2));
TST.rd=norm(TST.RES,'fro')/INV.normphi; %NOTE: unfiltered rd
TST.wrd = sum(rms(TST.RES) ./ (0.0010 + rms(INV.PHIREF))); % weighted rd

if strcmp(opt.type,'dep')
	TST.reg=norm(INV.REGOP*opt.tims)/1000;
else
	TST.reg=norm(INV.REGOPREP*opt.tims)/1000;
end
if INV.useWeighedRd
    TST.tres=sqrt(TST.wrd^2+(TST.reg*opt.mu)^2);
else
    TST.tres=sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
end

%% =======================================================================
function [REGOP,REGOPREP]=calcREGOP(GEOM)

if max(max(abs(GEOM.VER)))>1
    REGOP=surflapl(GEOM.VER/1000,GEOM.ITRI,1);
else
    REGOP=surflapl(GEOM.VER,GEOM.ITRI,1);
end
REGOPREP=REGOP;

%% =======================================================================
   
function S=getS(T,dep,rep,p,mode)

if mode == 4
    if size(T,2) < max(rep) + 50
        Tt=ones(length(dep),1)*([0:max(rep) + 50]);     
        S=getSmode(Tt,dep,rep,p,[],4);
        S=S(:,1:size(T,2));
    else
        S=getSmode(T,dep,rep,p,[],4);
    end
else
    S=getSmode(T,dep,rep,p,[],mode);
end