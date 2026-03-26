
function meas = inverse_pvd(varargin)%GEOM,inittimes,initcase,runcase,usetimes,leads,mode,mudep,murep)
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
INV.DOPLOT =0;
% inverse procedure parameters
INV.mudep           = 1.5e-4;
INV.murep           = 1.5e-4;
INV.muampl          = 1.5e-4;

INV.repOpt          = 'apd';
INV.MINRD           = 0.15;
INV.maxiter         = 25;
INV.lambopt         = 0.1;
INV.useAmpl         = 0;
INV.useNotch        = 0;
INV.useSurfLapl     = 1;
INV.useWeighedRd    = 0;
INV.experimentpvd   = 1;
INV.casedir         = 'invresults\';

INV.subname         = '';
INV.leads           = 1:65;
INV.mode            = 4;
INV.doWeight        = 0;
INV.doAmplitude     = 0;
INV.reg             = [];

% criteria for ending iterations
INV.lambdamax=500;
INV.stopcrit=2e-4;

%%
if length(varargin) < 3
    error('This routine needs at least three parameters');
else
    GEOM = varargin{1};
    INV.amplitude = ones(size(GEOM.VER,1),1);
    INV.initdep = varargin{2};
    INV.initrep = varargin{3};
    INV.leads = 1:length(GEOM.TVER);
    usetimes = size(GEOM.BSM,2) - GEOM.SPECS.onsetqrs;
    pp=4;
    while pp<=nargin
        if ischar(varargin{pp})
            key=lower(varargin{pp});
            switch key
                case 'mudep'
                    INV.mudep = varargin{pp+1};pp=pp+2;
                case 'murep'
                    INV.murep = varargin{pp+1};pp=pp+2;
                case 'repopt'
                    INV.repOpt = varargin{pp+1};pp=pp+2;
                case 'muampl'
                    INV.muampl = varargin{pp+1};pp=pp+2;
                case 'estimateampl'
                    INV.useAmpl = varargin{pp+1};pp=pp+2;
                case 'estimatedamplitude'
                    INV.amplitude = varargin{pp+1};pp=pp+2;
                case 'estimatenotch'
                    INV.useNotch = varargin{pp+1};pp=pp+2;
                case 'maxiter'
                    INV.maxiter = varargin{pp+1};pp=pp+2;
                case 'minrd'
                    INV.MINRD = varargin{pp+1};pp=pp+2;
                case 'casedir'
                    INV.casedir = varargin{pp+1};pp=pp+2;
                case 'logname'
                    INV.subname = varargin{pp+1};pp=pp+2;
                case 'usetimes'
                    INV.usetimes = varargin{pp+1};pp=pp+2;
                case 'leads'
                    INV.leads = varargin{pp+1};pp=pp+2;
                case 'weighed'
                    INV.doWeight = varargin{pp+1};pp=pp+2;
                case 'weighedrd'
                    INV.useWeighedRd  = varargin{pp+1};pp=pp+2;
                case 'mode'
                    INV.mode = varargin{pp+1};pp=pp+2;
                case 'reg'
                    INV.reg = varargin{pp+1};pp=pp+2;
                otherwise
                    error('unknown parameter');
            end
        end
    end
end

%%
if ~strcmp(INV.repOpt,'rep')
    INV.repOpt = 'apd';
end
INV.doSave = 0;
if ~exist(INV.casedir,'dir')
    mkdir(INV.casedir);
end
subname = [ '_' num2str(size(GEOM.LAY,1)-1)];
runcase = fullfile(INV.casedir,[ GEOM.subject GEOM.beat subname]); % oostep1
versie=1;


% PARAMETERS



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE INPUT DATA

% INV stucture
INV.lpass= 1;				% moving avarage lowpass filtering
INV.VER  = GEOM.VER;
INV.ITRI = GEOM.ITRI;
INV.ADJ  = GEOM.ADJ;
INV.DIST = GEOM.DIST;
INV.neigh = GEOM.neigh;

if strcmp(GEOM.type,'_atria') % no baseline correction for the atria
    
    if INV.mode == 1            
        INV.BSM =  lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetP:GEOM.SPECS.endtwave), INV.lpass);
        INV.BSM =  baselinecor(INV.BSM, GEOM.SPECS.onsetP, GEOM.SPECS.endP );
        INV.usetimes=min(GEOM.SPECS.endP-GEOM.SPECS.onsetP,size(INV.BSM,2));
    else        
        INV.BSM =  baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetP:GEOM.SPECS.endtwave), INV.lpass));
        INV.usetimes=min(usetimes,size(INV.BSM,2));
    end
    INV.PHIREF = INV.BSM(:,1:INV.usetimes);
else
    INV.BSM =  baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave), INV.lpass));
%     if INV.mode == 1
%         INV.usetimes = min(GEOM.SPECS.qrsduration,size(INV.BSM,2));
%     else
        INV.usetimes = min(usetimes,size(INV.BSM,2));
%     end
    INV.PHIREF = INV.BSM(:,1:INV.usetimes);
end

INV.normphi = norm(INV.PHIREF,'fro');
INV.T = ones(size(GEOM.AMA,2),1)*(0:INV.usetimes-1);

% compute regularization operator, e.g.surface Laplacian
if isempty(INV.reg)
    [INV.REGOP,INV.REGOPREP]=calcREGOP(GEOM,INV.useSurfLapl);
else
    INV.REGOP = INV.reg;
    INV.REGOPREP = INV.REGOP;
end
%precompute large matrixes
INV.AMA = GEOM.AMA;
INV.ATA = INV.AMA'*INV.AMA;

INV.weight=[];
if INV.doWeight
    bsm = GEOM.BSM(:,GEOM.SPECS.onsetqrs + GEOM.SPECS.time_Jpoint:end);
    INV.noiseLevel = mean(rms(bsm - lowpassma(bsm,21)));
    rrms=rms(INV.PHIREF);
    INV.weight = INV.noiseLevel./ rrms;
    INV.weight(rrms==0) = 0;
    INV.weight(500:end)=min(INV.weight(INV.weight>0));
    a = INV.weight(1:GEOM.SPECS.qrsduration) ./ (sum(INV.weight(1:GEOM.SPECS.qrsduration)));
    b = INV.weight(GEOM.SPECS.qrsduration + 1:end) ./ (sum(INV.weight(GEOM.SPECS.qrsduration+1:end)));
    INV.weight = [a b];
    INV.weight = INV.weight./mean(INV.weight(INV.weight > INV.weight(end)));
    INV.weight = 1+INV.weight./max(INV.weight);
    INV.weight = lowpassma(INV.weight,10);
    
    %     INV.weigth=INV.weight(:,1:size(INV.PHIREF,2))
    INV.normphiWeighed = norm(bsxfun(@times,INV.PHIREF,INV.weight),'fro'); % oostep1: 20121015 much faster
    %     weight = rms(INV.PHIREF);
    %     weight = weight/ max(weight);
    %     weight( weight < max(weight)/1.5) = max(weight)/1.5;
    %     INV.PHIREF = INV.PHIREF .* (weight * ones(1,size(INV.PHIREF,2)));
    %     INV.AMA = INV.AMA .*(weight * ones(1,size(INV.AMA,2)));
end

INV.ROTRO = INV.REGOP'*INV.REGOP;
INV.ROTROREP = INV.REGOPREP'*INV.REGOPREP;

INV.SPECS = GEOM.SPECS;
%% read initial estimates;

% regularization parameters; must be used to tune the final result:
% small values may produce spatially 'wild' solutions: use trial and error

% depolarization
OPT.DEP.mu      = INV.mudep; %1.4e-4
OPT.DEP.lambopt = INV.lambopt;
OPT.DEP.tims    = INV.initdep;
OPT.DEP.type    = 'dep';

% repolarization
OPT.REP.tims    = INV.initrep;
OPT.REP.lambopt = INV.lambopt;
OPT.REP.type    = INV.repOpt;
OPT.REP.mu      = INV.murep;

% Jpoint
OPT.NOT.pot=zeros(size(GEOM.VER,1),1);
if INV.useNotch% future when also an notch is included in the source model
    OPT.NOT.pot      = 1 - sqrAmpl(GEOM,GEOM.SPECS.time_Jpoint + 10 : GEOM.SPECS.time_Jpoint+40);
    OPT.AMP.usetimes = GEOM.SPECS.time_Jpoint:GEOM.SPECS.time_Jpoint + 60;
    OPT.NOT.lambopt  = INV.lambopt;
    OPT.NOT.mu       = INV.muampl;
    OPT.NOT.type     = 'notch';
end

% AP amplitude
OPT.AMP.pot = INV.amplitude;
if INV.useAmpl
    OPT.AMP.pot      = sqrAmpl(GEOM,GEOM.SPECS.time_Jpoint +0 : GEOM.SPECS.time_Jpoint + 20);
    OPT.AMP.usetimes = GEOM.SPECS.time_Jpoint:GEOM.SPECS.time_Jpoint + 20;
    OPT.AMP.lambopt  = INV.lambopt;
    OPT.AMP.mu       = INV.muampl;
    OPT.AMP.type     ='amplitude';
    doAmplitude      = 1;
end
%%
OPTSTART=OPT.DEP;
TST=gettres_v(INV,OPTSTART,OPT.REP,OPT.NOT,OPT.AMP);
%% =======================================================================
% DOCUMENT CASE
if INV.doSave
%     saveasci([runcase '.srcinit'],[initdep initrep]);
    logfile=[runcase '.log'];
    outfil= [runcase '.src'];
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
end
RESNOW=[];
%%
% PRELUDE TO ITERATIVE OMPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute reg and res and tres for initial estimate
if INV.doSave
    fprintf(fh,'%s\n','iter   min    mean     max     STD    reg     rd     tresd');
    temp=[0 min(OPT.DEP.tims) diff(range(OPT.DEP.tims)) max(OPT.DEP.tims) std(OPT.DEP.tims) TST.reg TST.rd TST.tres];
    sprintf(   '%3d %7.1f %7.1f %7.1f  %6.1f %6.1f %7.4f %7.4f',temp(1:8));
    fprintf(fh,'%3d %7.1f %7.1f %7.1f  %6.1f %6.1f %7.4f %7.4f\n',temp(1:8));

    fprintf(fh,'%s\n','iter   min   range  max    std  minAPD rangeAPD maxAPD stdAPD  reg    rd     tresd');
end
disp('iter   min   range  max    STD  minAPD rangeAPD maxAPD stdAPD  reg    rd     tresd');
temp=[0 min(OPT.DEP.tims) diff(range((OPT.DEP.tims))) max(OPT.DEP.tims) std(OPT.DEP.tims)...
    min(OPT.REP.tims-OPT.DEP.tims) diff(range(OPT.REP.tims-OPT.DEP.tims)) ...
    max(OPT.REP.tims-OPT.DEP.tims) std(OPT.REP.tims-OPT.DEP.tims)...
    TST.reg TST.rd TST.tres];
txt=sprintf(   '%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f Start time: %s',temp,datestr(now,'HH,MM.SS'));
disp(txt);
if INV.doSave
    fprintf(fh,'%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f',temp);
    fprintf(fh,'time run: %s\n',datestr(now,'HH,MM.SS')) ;
end
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
if INV.mode==-4, loopstart=1;nloop=2; end			% dep & rep ECGSIM
if INV.mode==5, loopstart=1;nloop=2; end			% dep & notch & rep
if INV.mode==6, loopstart=1;nloop=1; repscore=0;end % dep & notch
if INV.mode==7, loopstart=7;nloop=7; end			% AP amplitude only
totscore=0;
TST.rd=1;
while iter < 10 ||...
     iter < INV.maxiter && (depscore || repscore) && ...
     TST.rd > INV.MINRD
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
        if INV.doSave
            fprintf(fh,'%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f',temp);
            fprintf(fh,'run %s time %s\n',opt.type,datestr(datenum(clock)-datenum(startTime),'HH,MM.SS'));
        end
        RESNOW(ik,:)=[temp norm(TST.RES)];
        ik=ik+1;
    end
end
%%%% end outer loop; (loop1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TST=gettres_v(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);
meas.depfinal=OPT.DEP.tims;
meas.repfinal=OPT.REP.tims;
if INV.useAmpl
    meas.amplfinal=OPT.AMP.pot;
end
if INV.useNotch
    meas.notchPotfinal=OPT.NOT.pot;
end
meas.rdfinal=TST.rd;
COR=corrcoef(TST.PHIA,INV.PHIREF);
meas.corfinal=COR(2,1);
meas.iterfinal=iter;
meas.log=RESNOW;

if INV.doSave
    saveasci(outfil,[meas.depfinal meas.repfinal]);
    fclose(fh);
end

S=getS(INV.T,OPT.DEP.tims,OPT.REP.tims,INV.SPECS,OPT.NOT.pot,OPT.AMP.pot,INV.mode,INV.neigh);
PHIA=lowpassma(INV.AMA * S, INV.lpass);
TST.RES=INV.PHIREF - TST.PHIA(1:size(INV.PHIREF,1),1:size(INV.PHIREF,2));
if INV.DOPLOT
    figure(3);clf
    sigplot(INV.PHIREF,'',GEOM.LAY,1.5,'b',1,0);
    hold on
    sigplot(PHIA(1:size(INV.PHIREF,1),1:size(INV.PHIREF,2)),'',GEOM.LAY,1.5,'r',1,0);
end    

%% =======================================================================
function [score,amplopt,TST] = optimizeAmpl(INV,amplopt,depopt,repopt,notchopt)

% prepare compution Marquardt step aimed at improving tims;
% compute Sprime and Sprime*Sprime' based on previous iteration
TST=gettres_pot(INV,amplopt,depopt,repopt,notchopt);
starttres=TST.tres;
score=0;
lamb=amplopt.lambopt;

Splus=getS(INV.T,depopt.tims,repopt.tims,INV.SPECS,notchopt.pot,amplopt.pot+.01,4,INV.neigh);
Smin =getS(INV.T,depopt.tims,repopt.tims,INV.SPECS,notchopt.pot,amplopt.pot-.01,4,INV.neigh);
Sprime=(Splus(:,amplopt.usetimes)-Smin(:,amplopt.usetimes))/0.02;
SST=Sprime*Sprime';
% compute GTGM, with mu setting the weight of the regularization
% (smoothing based on the surface Laplacian)
GTGM=INV.ATA.*SST+amplopt.mu^2*INV.ROTRO;
% compute gtres
gtres=sum(INV.AMA.*(TST.RES*Sprime'))';
righths=gtres-amplopt.mu^2*INV.ROTRO*amplopt.pot;

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
TST=gettres_pot(INV,amplopt,depopt,repopt,notchopt);

% TST=gettres_v(INV,depopt,repopt,notchopt,amplopt);
%% =======================================================================
function [score,notchopt,TST]=optimizeNotch(INV,notchopt,depopt,repopt,amplopt)

% prepare compution Marquardt step aimed at improving tims;
% compute Sprime and Sprime*Sprime' based on previous iteration
TST=gettres_pot(INV,notchopt,depopt,repopt,amplopt);
starttres=TST.tres;
score=0;
lamb=notchopt.lambopt;

Splus=getS(INV.T,depopt.tims,repopt.tims,INV.SPECS,notchopt.pot+0.01,amplopt.pot,4,INV.neigh);
Smin =getS(INV.T,depopt.tims,repopt.tims,INV.SPECS,notchopt.pot-0.01,amplopt.pot,4,INV.neigh);
Sprime=(Splus(:,notchopt.usetimes)-Smin(:,notchopt.usetimes)) / 0.02;
SST=Sprime*Sprime';

% compute GTGM, with mu setting the weight of the regularization
% (smoothing based on the surface Laplacian)
M = INV.ATA.*SST;
GTGM = M'*M + + notchopt.mu^2*INV.ROTRO;
% GTGM = M + muREG; old
gtres = sum(INV.AMA.*(TST.RES*Sprime'))';

righths= M' * gtres - muREG * startopt.tims;
righths = gtres - notchopt.mu^2*INV.ROTRO*notchopt.pot;
	    

testopt=startopt;
% compute next estimate; test residual;
% increase lambda it until starttres < TST.tres

% compute gtres
righths=gtres - notchopt.mu^2*INV.ROTRO*notchopt.pot;

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
TST=gettres_pot(INV,notchopt,depopt,repopt,amplopt);

% TST=gettres_v(INV,depopt,repopt,notchopt,amplopt);

%% =======================================================================
function [score,opt,TST]=optimizeDepRep(INV,opt,keepopt,notchopt,amplopt)

% prepare compution Marquardt step aimed at improving tims;
% compute Sprime and Sprime*Sprime' based on previous iteration

TST=gettres_v(INV,opt,keepopt,notchopt,amplopt);
starttres=TST.tres;
score=0;
lamb=opt.lambopt;
delta = 1;
if strcmp(opt.type,'dep')
    delta = 1;
    SPECS  = INV.SPECS;
    SPECS.depSlope = 1;
    Splus=getS(INV.T,opt.tims+delta,keepopt.tims,SPECS,notchopt.pot,amplopt.pot,INV.mode,INV.neigh);
    Smin =getS(INV.T,opt.tims-delta,keepopt.tims,INV.SPECS,notchopt.pot,amplopt.pot,INV.mode,INV.neigh);
    startopt=opt;
elseif strcmp(opt.type,'rep')
    delta = 1;
    Splus=getS(INV.T,keepopt.tims,opt.tims+delta,INV.SPECS,notchopt.pot,amplopt.pot,INV.mode,INV.neigh);
    Smin =getS(INV.T,keepopt.tims,opt.tims-delta,INV.SPECS,notchopt.pot,amplopt.pot,INV.mode,INV.neigh);
    startopt=opt;
else %if strcmp(opt.type,'apd')
    delta = 1;
    Splus=getS(INV.T,keepopt.tims,opt.tims+delta,INV.SPECS,notchopt.pot,amplopt.pot,INV.mode,INV.neigh);
    Smin =getS(INV.T,keepopt.tims,opt.tims-delta,INV.SPECS,notchopt.pot,amplopt.pot,INV.mode,INV.neigh);
    startopt=opt;
    startopt.tims=startopt.tims-keepopt.tims;
end

Sprime=(Splus-Smin)/(2*delta);
SST=Sprime*Sprime';
% compute GTGM, with mu setting the weight of the regularization
% (smoothing based on the surface Laplacian)
if strcmp(startopt.type,'dep')
    % experiment pvd
    if INV.experimentpvd
        dep = opt.tims;
        dep =abs(dep- mean(dep));
        dep = dep /max(dep);
        mu = (1-(dep.^2));
        %     mu(mu < 0.5) = mu(mu < 0.5)/5;
        mu = (opt.mu^2).* mu;
        %     mu = (opt.mu^2).*((dep.^2));
        muREG = bsxfun(@times,mu,INV.ROTRO);
        M = INV.ATA.*SST;
        
    else
        muREG = bsxfun(@times,opt.mu^2,INV.ROTRO);
        M = INV.ATA.*SST;
        
    end
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
    if lamb > INV.lambdamax,  
        break;
    end
    GTGM_L = GTGM + lamb^2*eye(size(INV.AMA,2));
%     lastwarn('','');
%     invGTGML = inv(GTGML);
%     if ~isempty(lastwarn)
%         disp('pinv used');
%         invGTGML =pinv(GTGM_L);
%     end
%     deltau = invGTGML * righths;
        
    deltau=GTGM_L \ righths;  % INVERSE; total computing time  0.87 times preceding

    newtime=startopt.tims+deltau;
%     if strcmp(startopt.type,'dep')
%         newtime(newtime<-5)=-5;
%     elseif strcmp(startopt.type,'rep')
%         newtime(newtime<max(keepopt.tims)+50)=max(keepopt.tims)+50;
%     elseif strcmp(startopt.type,'apd')
%         newtime(newtime<50)=50;
%     end
    
    testopt.tims=newtime;
    if strcmp(startopt.type,'apd')
        testopt.tims = testopt.tims + keepopt.tims;
        TST=gettres_v(INV,testopt,keepopt,notchopt,amplopt);
    else
        TST=gettres_v(INV,testopt,keepopt,notchopt,amplopt);
    end
    if TST.tres < starttres
        if (starttres-TST.tres)/starttres >= INV.stopcrit
            opt.tims=testopt.tims;
            score=1;
        end
        break;
    end
end  % end of inner loop
opt.lambopt=lamb/4;

%% =======================================================================
function TST = gettres_pot(INV,opt,depopt,repopt,otheropt)

if strcmp(opt.type,'notch')
    S=getS(INV.T,depopt.tims,repopt.tims,INV.SPECS,opt.pot,otheropt.pot,INV.mode,INV.neigh);
else %type == ' ampl',
    S=getS(INV.T,depopt.tims,repopt.tims,INV.SPECS,otheropt.pot,opt.pot,INV.mode,INV.neigh);
end
TST.PHIA=lowpassma(INV.AMA*S(:,opt.usetimes), INV.lpass);
TST.RES=INV.PHIREF(:,opt.usetimes) - TST.PHIA;
TST.rd=norm(TST.RES,'fro')/norm(INV.PHIREF(:,opt.usetimes),'fro'); %NOTE: unfiltered rd
disp('remove 1000?');
% TST.reg=norm(INV.REGOP*opt.pot)/1000;
TST.reg=norm(INV.REGOP*opt.pot);
TST.tres=sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
TST.S=S;

%% =======================================================================
function TST=gettres_v(INV,opt,keepopt,notchopt,amplopt)
% global temp
if strcmp(opt.type,'dep')
    TST.S=getS(INV.T,opt.tims,keepopt.tims,INV.SPECS,notchopt.pot,amplopt.pot,INV.mode,INV.neigh);
else %pol=='rep',
    TST.S=getS(INV.T,keepopt.tims,opt.tims,INV.SPECS,notchopt.pot,amplopt.pot,INV.mode,INV.neigh);
end


TST.PHIA=lowpassma(INV.AMA*TST.S, INV.lpass);
TST.RES = INV.PHIREF - TST.PHIA(1:size(INV.PHIREF,1),1:size(INV.PHIREF,2));


if INV.mode == 1 %|| strcmp(opt.type,'dep')
    maxt = ceil(min([max(opt.tims),size(INV.PHIREF,2)]) );
    normphi = norm(INV.PHIREF(:,1:maxt),'fro');
    TST.rd=norm(TST.RES(:,1:min(size(TST.RES,2),ceil(max(opt.tims)))),'fro')/normphi; %NOTE: unfiltered rd
else
    normphi  = INV.normphi;
    TST.rd=norm(TST.RES,'fro')/INV.normphi; %NOTE: unfiltered rd
end
% scale = 1./ (0.10 + rms(INV.PHIREF));

% TST.wrd = norm(bsxfun(@times,TST.RES,scale),'fro') ./ norm(bsxfun(@times,INV.PHIREF,scale),'fro');
% TST.wrd = norm(bsxfun(@times,TST.RES,scale),'fro') ./ normphi;
% 
% 
% TST.wrd = sum(rms(TST.RES) ./ (0.010 + rms(INV.PHIREF))); % weighted rd



if strcmp(opt.type,'dep')
%     velo=calcActVelo(opt.tims,INV.ADJ,INV.DIST,find(opt.tims==min(opt.tims)));
%     velo=calcActVelo(opt.tims,INV.DIST,30);
%     TST.reg=norm(INV.REGOP*velo);
%  L=surflapl(INV.VER,INV.ITRI,0);area=nodearea(INV.VER,INV.ITRI);a=sqrt(sum((L*opt.tims).^2.*(area)))
    TST.reg=norm(INV.REGOP * opt.tims);
else
    TST.reg=norm(INV.REGOPREP*opt.tims);
end
if INV.useWeighedRd
    if strcmp(opt.type,'dep')
%         disp(num2str([TST.wrd (TST.reg*opt.mu)]))
        TST.tres=sqrt(TST.wrd^2 + (TST.reg*opt.mu)^2);
    else
        TST.tres=sqrt(TST.wrd^2 + (TST.reg*opt.mu)^2);
    end
else
    if strcmp(opt.type,'dep')
%         disp(num2str([TST.rd (TST.reg*opt.mu) TST.reg]))
        
        TST.tres=sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
    else
        TST.tres=sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
    end
end
%% =======================================================================
function [REGOP,REGOPREP]=calcREGOP(GEOM,useSurfLapl)

% if isfield(GEOM,'LAPL')
%     REGOP2 = 2.0 * GEOM.LAPL;
%     REGOPREP2 = REGOP;
% end
if useSurfLapl==1
    if max(max(abs(GEOM.VER)))>1
%         REGOP=surflapl(GEOM.VER/1000,GEOM.ITRI,1);
        REGOP=surflapl(GEOM.VER,GEOM.ITRI,1);
    else
        REGOP=surflapl(GEOM.VER,GEOM.ITRI,1);
    end
    %     if isfield(GEOM,'LAPL')
    %         REGOP = REGOP - GEOM.LAPL;
    %     end
    REGOPREP=REGOP;
    %
    %     REGOPREP=calcAnisADJ(GEOM,0.25);
    % 	REGOPREP(GEOM.ADJsurf>0)=0;
    % 	for i=1:length(REGOPREP)
    % 		a=[(1:length(REGOPREP))' REGOPREP(:,i)];a=sortrows(a,2);
    % 		a(a(:,2)==0,:)=[];a=a(8:end,:);
    % 		REGOPREP(a(:,1),i)=0;		REGOPREP(i,a(:,1))=0;
    % 	end
    % 	REGOPREP(GEOM.ADJsurf>0)=GEOM.ADJsurf(GEOM.ADJsurf>0);
    % % 	REGOPREP=graphdist(REGOPREP);
    % 	if max(max(abs(GEOM.VER)))>1
    % 		REGOPREP=REGOPREP/1000;
    % 	end
    % 	REGOPREP=REGOPREP.^(-2);
    % 	REGOPREP(REGOPREP==Inf)=0;
    % 	scal=sum(REGOPREP,2);
    % 	scal=1./scal;
    % 	REGOPREP=diag(scal)*REGOPREP;
    % 	REGOPREP=(REGOPREP-eye(size(REGOPREP)));% scales the result to be of the same order of SL
    % 	fact=[];
    % 	for i=1:length(REGOPREP)
    % 		fact(i)=REGOP(i,i)/REGOPREP(i,i);
    % 	end
    % 	REGOPREP=REGOPREP*mean(fact);
    %     % 59.7
    % 	REGOPREP=597*(REGOPREP-eye(size(REGOPREP)));% scales the result to be of the same order of SL
    % REGOP(1,REGOP(1,:)~=0)
    % REGOPREP(1,REGOPREP(1,:)~=0)
    % REGOPREP=REGOP;
else %if useSurfLapl==2
    disp('inverse 1.distance^2 regularization');
    [a,ORDER] = graphdist(GEOM.ITRI);
    REGOP = GEOM.DIST2W;
    %     REGOP=GEOM.ANISDIST;
    %     ADJ=GEOM.ADJ; ADJ(ADJ>50)=0;
    %     REGOP=calcAnisADJ(GEOM.VER,GEOM.ITRI,ADJ,GEOM.RVER,GEOM.LVER,2.5);
    %     REGOP=calcAnisADJ(GEOM,2.5);
    if max(max(abs(GEOM.VER)))>1
        REGOP=3*REGOP/1000;
    end
    REGOP=REGOP.^(-1);
    REGOPREP = REGOP;
    if max(max(abs(GEOM.VER)))>1
%         B=surflapl(GEOM.VER/1000,GEOM.ITRI,1);
        B=surflapl(GEOM.VER,GEOM.ITRI,1);
    else
        B=surflapl(GEOM.VER,GEOM.ITRI,1);
    end
    
    REGOP(GEOM.DIST>10)=0;
    REGOP(ORDER<=4)=0;
    B2=diag(-0.3*sum(REGOP));
    REGOP(B~=0)=B(B~=0);
    REGOP =  REGOP + B2;
    %     REGOPREP(GEOM.DIST>20)=0;
    %     REGOPREP(ORDER<=1)=0;
    %     REGOPREP(B~=0)=B(B~=0);
    REGOPREP=B;
    
    % 	REGOP(REGOP==Inf)=0;
    % 	scal=sum(REGOP,2);
    % 	scal=1./scal;
    % 	REGOP=diag(scal)*REGOP;
    %     % 59.7
    % 	REGOP=70*(REGOP-eye(size(REGOP)));% scales the result to be of the same order of SL
end

%% =======================================================================

function S=getS(T,dep,rep,SPECS,notchpot,scaleAmpl,mode,neigh)

if mode == 4 || mode == 5
    SPECS.scaleAmpl = scaleAmpl;
    if size(T,2) < max(rep) + 50
        Tt=ones(length(dep),1)*(0:max(rep) + 50);
        S=getSmode(Tt,dep,rep,SPECS,4,neigh);
        S=S(:,1:size(T,2));
    else
        S=getSmode(T,dep,rep,SPECS,4,neigh);
    end
else
    SPECS.scaleAmpl = notchpot;
    S=getSmode(T,dep,rep,SPECS,mode,neigh);
end
