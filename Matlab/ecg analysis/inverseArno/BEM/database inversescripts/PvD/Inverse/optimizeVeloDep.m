
function meas=optimizeVeloDep(varargin)%GEOM,inittimes,initcase,runcase,usetimes,leads,mode,mudep,murep)
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

MINRD=0.15;
maxiter=10; 
lambopt=0.1;  
useAmpl = 0;
useNotch =0;
casedir	='invresults\';
subname='';
leads = 1:64;
mode = 1;
doWeight = 0;
%%
if length(varargin) < 3
	error('This routine needs at least three parameters');
else
    GEOM = varargin{1};
%     initdep = varargin{2};
    foci = varargin{2};
    fociTimes = varargin{3};
    usetimes = size(GEOM.BSM,2)-GEOM.specs(2)+1;
	pp=4;
    while pp<=nargin
        if ischar(varargin{pp})
            key=lower(varargin{pp});
            switch key
                case 'mudep'
                    mudep = varargin{pp+1};pp=pp+2;
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
                otherwise
                    error('unknown parameter');
            end
        end
    end
end


if ~exist(casedir,'dir')
    mkdir(casedir);
end
runcase = [casedir GEOM.subject GEOM.beat subname];
versie=1;	


% PARAMETERS
lpass=5;				% moving avarage lowpass filtering
useSurfLapl=0;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE INPUT DATA

% INV stucture
INV.VER=GEOM.VER;
INV.ITRI=GEOM.ITRI;
INV.foci = foci;
INV.fociTimes = fociTimes;


INV.mode=mode;

if strcmp(GEOM.type,'_atria') % no baseline coroection for the atria
    INV.BSM=zeromean(GEOM.BSM(:,GEOM.specs(2):end));
    if INV.mode == 1
         INV.usetimes=min(GEOM.specs(3)-GEOM.specs(2)+1,size(INV.BSM,2));
    else
        INV.usetimes=min(usetimes,size(INV.BSM,2));
    end
    INV.PHIREF=zeromean(INV.BSM(leads,1:INV.usetimes));
else
    INV.BSM=baselinecor(zeromean(GEOM.BSM(:,GEOM.specs(2):end)));
    if INV.mode == 1
        INV.usetimes=min(GEOM.specs(3)-GEOM.specs(2)+1,size(INV.BSM,2));
    else
        INV.usetimes=min(usetimes,size(INV.BSM,2));
    end
    INV.PHIREF=baselinecor(zeromean(INV.BSM(leads,1:INV.usetimes)));
end

INV.ADJ   = GEOM.ADJ;
INV.ADJ2W = GEOM.ADJ2W;


INV.normphi=norm(INV.PHIREF,'fro');
INV.T=ones(size(GEOM.AMA,2),1)*(0:INV.usetimes-1);

% compute regularization operator, e.g.surface Laplacian
INV.REGOP=calcREGOP(GEOM,useSurfLapl);
%precompute large matrixes
INV.AMA=zeromean(GEOM.AMA(leads,:));
INV.ATA=INV.AMA'*INV.AMA;  

if doWeight
   weight = rms(INV.PHIREF')';
    weight = weight/ max(weight);
    weight( weight < max(weight)/1.5) = max(weight)/1.5;
    INV.PHIREF = INV.PHIREF .* (weight * ones(1,size(INV.PHIREF,2)));
    INV.AMA = INV.AMA .*(weight * ones(1,size(INV.AMA,2)));
end
INV.REG = INV.REGOP'*INV.REGOP;


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
% OPT.DEP.tims=initdep;
OPT.DEP.type='dep';
OPT.DEP.velocity = 0.8*ones(length(INV.ADJ2W),1);

%%
[TST,OPT.DEP]=gettres_velo(INV,OPT.DEP);
%% =======================================================================
% DOCUMENT CASE  
% saveasci([runcase 'srcinit'],[initdep initrep]);
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
fprintf(fh,'muvelo: %0.4g   \n',OPT.DEP.mu);
fprintf(fh,'%s\n','ifix(1) tfix(1) dampd   dampr usetimes   tblc   mode  maxiter');

RESNOW=[];
%%    
% PRELUDE TO ITERATIVE OMPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute reg and res and tres for initial estimate

fprintf(fh,'%s\n','iter   min    mean     max     STD    reg     rd     tresd');
temp=[0 min(OPT.DEP.tims) diff(range(OPT.DEP.tims)) max(OPT.DEP.tims) std(OPT.DEP.tims) TST.reg TST.rd TST.tres];
sprintf(   '%3d %7.1f %7.1f %7.1f  %6.1f %6.1f %7.4f %7.4f',temp(1:8));
fprintf(fh,'%3d %7.1f %7.1f %7.1f  %6.1f %6.1f %7.4f %7.4f\n',temp(1:8));

fprintf(fh,'%s\n','iter   min   range  max    std  minVelo rangeVelo maxVelo stdVelo  reg    rd     tresd');
txt=sprintf('%s\n','iter   min   range  max    STD  minVelo rangeVelo maxVelo stdVelo  reg    rd     tresd');
disp(txt);
temp=[0 min(OPT.DEP.tims) diff(range((OPT.DEP.tims))) max(OPT.DEP.tims) std(OPT.DEP.tims)...
        min(OPT.DEP.velocity) diff(range(OPT.DEP.velocity)) max(OPT.DEP.velocity) std(OPT.DEP.velocity)...
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
depscore=1;

TST.rd=1;
while iter<maxiter && depscore && TST.rd>MINRD 
	iter=iter+1;
	
    [depscore,OPT.DEP,TST]=optimizeDepVelo(INV,OPT.DEP);
				
    temp=[iter min(OPT.DEP.tims) diff(range(OPT.DEP.tims)) max(OPT.DEP.tims) std(OPT.DEP.tims)...
          min(OPT.DEP.velocity) diff(range(OPT.DEP.velocity)) max(OPT.DEP.velocity) std(OPT.DEP.velocity)...
          TST.reg TST.rd TST.tres];
    if score==0
        k=rem(ik,2)+1;
        temp(9:10)=RESNOW(ik-k,9:10);
        temp(11)=RESNOW(ik-1,11);
    end
    txt=sprintf(   '%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f run %s time %s',...
                temp,OPT.DEP.type,datestr(datenum(clock)-datenum(startTime),'HH,MM.SS'));
    disp(txt);
    fprintf(fh,'%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f',temp);
    fprintf(fh,'run %s time %s\n',OPT.DEP.type,datestr(datenum(clock)-datenum(startTime),'HH,MM.SS'));
	RESNOW(ik,:)=[temp norm(TST.RES)];
	ik=ik+1;
end
%%%% end outer loop; (loop1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[TST,OPT.DEP]=gettres_velo(INV,OPT.DEP);
meas.depfinal=OPT.DEP.tims;
meas.rdfinal=TST.rd;
COR=corrcoef(TST.PHIA,INV.PHIREF);
meas.corfinal=COR(2,1);
meas.iterfinal=iter;
meas.log=RESNOW;

saveasci(outfil,[meas.depfinal meas.repfinal]);
saveasci(outfilP,INV.pSn);
fclose(fh);


%% =======================================================================
function [score,opt,TST]=optimizeDepVelo(INV,opt)

% prepare compution Marquardt step aimed at improving tims; 
% compute Sprime and Sprime*Sprime' based on previous iteration

[TST,opt]=gettres_velo(INV,opt);
starttres=TST.tres;
score=0;
lamb=opt.lambopt; 


[Splus,depPlus]=getS(INV,opt.velocity + 0.05);
[Smin,depMin] =getS(INV,opt.velocity - 0.05);
dt= depPlus - depMin;

A = zeros(size(depPlus));
for i=1:length(depPlus)
    a=find( INV.ADJ2W(i,:) > 0);
    A(i) = max(-0.1 ./ dt(a));
end

startopt=opt;	

Sprime=(Splus-Smin)/2.;
SST=Sprime*Sprime';
% compute GTGM, with mu setting the weight of the regularization
% (smoothing based on the surface Laplacian)
GTGM=INV.ATA.*SST+opt.mu^2*INV.REG;
% compute gtres
gtres=sum(INV.AMA.*(TST.RES*Sprime'))';
righths=gtres-opt.mu^2*INV.REG*startopt.velocity;

testopt=startopt;
% compute next estimate; test residual;
% increase lambda it until starttres < TST.tres

while TST.tres >= starttres, % inner loop
	lamb=lamb*2;
	if lamb > INV.lambdamax,  break, end
	GTGML = GTGM + lamb^2*eye(size(INV.AMA,2));
    
	deltVelo = pinv(GTGML)*righths;
    deltVelo = deltVelo /50;
    
	newVelocity = startopt.velocity + deltVelo; 
    newVelocity(newVelocity < 0.2)=0.2;
	
	
	testopt.velocity=newVelocity;
    [TST,testopt]=gettres_velo(INV,testopt);		
	if TST.tres<starttres
		if (starttres-TST.tres)/starttres >= INV.stopcrit
			opt.velocity = testopt.velocity; 
            opt.tims = testopt.tims;
			score=1;
		end
		break;
	end
end  % end of inner loop
opt.lambopt=lamb/4;

%% =======================================================================
function [TST,opt]=gettres_velo(INV,opt)

global lpass
[TST.S,opt.tims]   = getS(INV,opt.velocity);
TST.PHIA = lowpassma(INV.AMA*TST.S,lpass);
TST.RES  = INV.PHIREF-TST.PHIA(1:size(INV.PHIREF,1),1:size(INV.PHIREF,2));
TST.rd   = norm(TST.RES,'fro')/INV.normphi; %NOTE: unfiltered rd
TST.reg = norm(INV.REGOP*opt.tims)/1000;
TST.tres = sqrt(TST.rd^2+(TST.reg*opt.mu)^2);

%% =======================================================================
function REGOP=calcREGOP(GEOM,useSurfLapl)

if useSurfLapl==1
    if max(max(abs(GEOM.VER)))>1
        REGOP=surflapl(GEOM.VER/1000,GEOM.ITRI,1);
    else
        REGOP=surflapl(GEOM.VER,GEOM.ITRI,1);
    end
else %if useSurfLapl==2
	disp('inverse 1.distance^2 regularization');
    LD =GEOM.ADJ2W;
    LD(LD > 0) = (LD( LD > 0 ).^-1.5);
    REGOP = LD - diag(sum(LD));
end

%% =======================================================================
   
function varargout = getS(INV,velo)

ADJ = INV.ADJ2W/0.6;
for i=1:length(ADJ)
    ADJ(i,:) = ADJ(i,:) * (0.6 / velo(i));
end
dep = focipaths(ADJ/2,INV.foci,INV.fociTimes);

TDEP = INV.T - dep * ones(1,size(INV.T,2));
S = 1 ./ ( 1 + exp(2 * TDEP));	


varargout{1} = diag(1./max(S,[],2))*S;


if nargout > 1
    varargout{2} = dep;    
end
