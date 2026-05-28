function meas = inverse_aj(varargin) 
%
% CHECK: diff(range((OPT.DEP.tims)))
%
% INPUT:
% varargin{1} = GEOM;
% varargin{2} = INV.initdep;
% varargin{3} = INV.initrep;
% 
% OPTIONAL INPUT:
% 'mudep',          INV.mudep
% 'murep',          INV.murep
% 'repopt',         INV.repOpt
% 'muampl',         INV.muampl
% 'estimateampl',   INV.useAmpl
% 'maxiter',        INV.maxiter
% 'minrd',          INV.MINRD
% 'casedir',        INV.casedir
% 'logname',        INV.subname
% 'usetimes',       INV.usetimes
% 'leads',          INV.leads
% 'weighedrd',      INV.useWeighedRd
% 'mode',           INV.mode
% 'reg',            INV.reg
% 'blmode',         INV.blmode
% 'lpass',          INV.lpass
%
% inverse determination of timing of Equivalent Double Layer source EDL;
% optimization based on basic Marquardt type of solving the non-linear
% estimation problem alternate iterations between depolarization and
% repolarization.
%
% Peter van Dam; 2010 november.
% All rights reserved Peacs, Arnhem  the Netherlands
%
% if INV.mode == 1, loopstart = 1; nloop = 1; repscore = 0; end % dep only
% if INV.mode == 2, loopstart = 2; nloop = 2; depscore = 0; end % rep only
% if INV.mode == 3, loopstart = 3; nloop = 3; end % notch
% if INV.mode == 4, loopstart = 1; nloop = 2; end % dep rep
% if INV.mode == 5, loopstart = 1; nloop = 3; end % dep rep & notch
% if INV.mode == 6, loopstart = 6; nloop = 6; end % AP amplitude only
%
% Copy of inverse_po.m 9-jan-2015: version 1.0 AJ

%% inverse procedure parameters

% STANDARD PARAMETERS
INV.mudep           = 1.5e-4;           % depolarization
INV.murep           = 1.5e-4;           % repolarization
INV.muampl          = 1.5e-4;           %
INV.repOpt          = 'apd';            %
INV.MINRD           = 0.15;             % minimum RD
INV.maxiter         = 25;               % maximum number of iterations 
INV.lambopt         = 0.001;            % 
INV.useAmpl         = 0;                %
INV.useNotch        = 0;                %
INV.useSurfLapl     = 1;                %
INV.useWeighedRd    = 1;                %
INV.experimentpvd   = 0;                %
INV.casedir         = 'invresults\';    %
INV.doSave          = 0;                %
INV.subname         = '';               %
INV.leads           = 1:65;             %
INV.mode            = 4;                %
INV.doWeight        = 0;                %
INV.doAmplitude     = 0;                %
INV.reg             = [];               %
INV.blmode          = 1;                % [1] blcorrection qrsonset - T, [2] (Vstim-10) - T.
INV.lpass           = 5;                % moving avarage lowpass filtering
INV.lambdamax       = 500;              % criteria for ending iterations
INV.stopcrit        = 2e-4;             % criteria for ending iterations

% ADJUST PARAMETERS WITH INPUT:
if length(varargin) < 3,
    error('This routine needs at least three parameters');
else
    GEOM            = varargin{1};
    INV.amplitude   = ones(size(GEOM.VER,1),1);
    INV.initdep     = varargin{2};
    INV.initrep     = varargin{3};
    INV.leads       = 1:length(GEOM.TVER);
    usetimes        = size(GEOM.BSM,2) - GEOM.SPECS.onsetqrs;
    pp              = 4;
    while pp <= nargin
        if ischar(varargin{pp})
            key = lower(varargin{pp});
            switch key
                case 'mudep', INV.mudep              = varargin{pp+1}; pp = pp+2;
                case 'murep', INV.murep              = varargin{pp+1}; pp = pp+2;
                case 'repopt', INV.repOpt            = varargin{pp+1}; pp = pp+2;
                case 'muampl', INV.muampl            = varargin{pp+1}; pp = pp+2;
                case 'estimateampl', INV.useAmpl     = varargin{pp+1}; pp = pp+2;
                case 'maxiter', INV.maxiter          = varargin{pp+1}; pp = pp+2;
                case 'minrd', INV.MINRD              = varargin{pp+1}; pp = pp+2;
                case 'casedir', INV.casedir          = varargin{pp+1}; pp = pp+2;
                case 'logname', INV.subname          = varargin{pp+1}; pp = pp+2;
                case 'usetimes', INV.usetimes        = varargin{pp+1}; pp = pp+2;
                case 'leads', INV.leads              = varargin{pp+1}; pp = pp+2;
                case 'weighedrd', INV.useWeighedRd   = varargin{pp+1}; pp = pp+2;
                case 'mode', INV.mode                = varargin{pp+1}; pp = pp+2;
                case 'reg', INV.reg                  = varargin{pp+1}; pp = pp+2;
                case 'blmode', INV.blmode            = varargin{pp+1}; pp = pp+2;
                case 'lpass', INV.lpass              = varargin{pp+1}; pp = pp+2;
                otherwise, error('unknown parameter');
            end
        end
    end
end

if ~strcmp(INV.repOpt,'rep'), INV.repOpt = 'apd'; end
if ~exist(INV.casedir,'dir'), mkdir(INV.casedir); end   % check if output directory is present, else create

subname = [ '_' num2str(size(GEOM.LAY,1)-1)];
runcase = fullfile(INV.casedir,[GEOM.subject GEOM.beat subname]); % name of file
versie  = 1;

%% PREPARE INPUT DATA			
INV.VER     = GEOM.VER;
INV.ITRI    = GEOM.ITRI;
INV.ADJ     = GEOM.ADJ;
INV.DIST    = GEOM.DIST;
INV.neigh   = GEOM.neigh;

% BASELINE CORRECTION FOR DATA [BSM]:
if strcmp(GEOM.type,'_atria'), % no baseline correction for the atria
    if INV.blmode ~= 1, warning('Only baselinemode 1 supported for atria'); end 
    INV.BSM =  baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetP:GEOM.SPECS.endtwave),INV.lpass));
    if INV.mode == 1, INV.usetimes = min(GEOM.SPECS.onsetqrs,size(INV.BSM,2)); else INV.usetimes = min(usetimes,size(INV.BSM,2)); end
    INV.PHIREF = INV.BSM(:,1:INV.usetimes);
else
    if INV.blmode == 1, % order of baselinecor and lowpassma
        INV.BSM = baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),INV.lpass));
    else
        bsmt    = baselinecor(GEOM.BSM,GEOM.SPECS.time_Vstim-10,GEOM.SPECS.endtwave);
        INV.BSM = lowpassma(bsmt(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),5);
    end
    
    if INV.mode == 1, INV.usetimes = min(GEOM.SPECS.qrsduration,size(INV.BSM,2)); else INV.usetimes = min(usetimes,size(INV.BSM,2)); end
    
    INV.PHIREF          = INV.BSM(:,1:INV.usetimes);
end

INV.normphi = norm(INV.PHIREF,'fro');                       % Frobenius norm of BSM (during QRS)
INV.T       = ones(size(GEOM.AMA,2),1)*(0:INV.usetimes-1);  % matrix of nr of vertices on ventricles mesh x duration BSM (during QRS)

% Compute regularization operator, e.g.surface Laplacian
if isempty(INV.reg), 
    [INV.REGOP,INV.REGOPREP]    = calcREGOP(GEOM,INV.useSurfLapl);
else
    INV.REGOP                   = INV.reg;
    INV.REGOPREP                = INV.REGOP;
end

%precompute large matrixes
INV.AMA         = GEOM.AMA;             % Transfer (A) Matrix
INV.ATA         = INV.AMA'*INV.AMA;     % transpose-A-matrix * A-matrix

INV.ROTRO       = INV.REGOP'*INV.REGOP;
INV.ROTROREP    = INV.REGOPREP'*INV.REGOPREP;
INV.SPECS       = GEOM.SPECS;

%% read initial estimates;
% regularization parameters must be used to tune the final result:
% small values may produce spatially 'wild' solutions: use trial and error

% depolarization
OPT.DEP.mu      = INV.mudep;
OPT.DEP.lambopt = INV.lambopt;
OPT.DEP.tims    = INV.initdep;
OPT.DEP.type    = 'dep';

% repolarization
OPT.REP.mu      = INV.murep;
OPT.REP.lambopt = INV.lambopt;
OPT.REP.tims    = INV.initrep;
OPT.REP.type    = INV.repOpt;

% Jpoint
OPT.NOT.pot     = zeros(size(GEOM.VER,1),1);
if INV.useNotch, % future when also a notch is included in the source model
    OPT.NOT.pot      = 1-sqrAmpl(GEOM,GEOM.SPECS.time_Jpoint+10:GEOM.SPECS.time_Jpoint+40);
    OPT.AMP.usetimes = GEOM.SPECS.time_Jpoint:GEOM.SPECS.time_Jpoint+60;
    OPT.NOT.mu       = INV.muampl;
    OPT.NOT.lambopt  = INV.lambopt;
    OPT.NOT.type     = 'notch';
end

% AP amplitude
OPT.AMP.pot = INV.amplitude;
if INV.useAmpl,
    OPT.AMP.pot      = sqrAmpl(GEOM,GEOM.SPECS.time_Jpoint+10:GEOM.SPECS.time_Jpoint+40);
    OPT.AMP.usetimes = GEOM.SPECS.time_Jpoint:GEOM.SPECS.time_Jpoint+60;
    OPT.AMP.lambopt  = INV.lambopt;
    OPT.AMP.mu       = INV.muampl;
    OPT.AMP.type     ='amplitude';
    doAmplitude      = 1;
end

TST         = gettres_v(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);

%% =======================================================================
% DOCUMENT CASE
if INV.doSave,
    logfile = [runcase '.log'];
    outfil  = [runcase '.src'];
    fh      = fopen(logfile,'w');
    fprintf(fh,'%s\n',['file: ' logfile]);
    fprintf(fh,'%s\n',['date: ' datestr(clock)]);
    fprintf(fh,'%s\n',['prog: invgeom (' num2str(versie) ')' ]);
    fprintf(fh,'%s\n',['initvals:   ' [runcase 'srcinit']]);
    fprintf(fh,'%s\n',['subject:   ' GEOM.subject]);
    fprintf(fh,'%s\n',['used anisotropyRatio:    ' GEOM.anisotropyRatio]);
    
    if useSurfLapl == 1, 
        fprintf(fh,'%s\n','SL regularization');
    else
        fprintf(fh,'%s\n','inverse 1.distance^2 regularization ');
    end
    
    fprintf(fh,'mudep: %0.4g    murep: %0.4g \n',OPT.DEP.mu,OPT.REP.mu);
    fprintf(fh,'%s\n','ifix(1) tfix(1) dampd   dampr usetimes   tblc   mode  maxiter');
end

RESNOW = [];

%% PRELUDE TO ITERATIVE OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute reg and res and tres for initial estimate
if INV.doSave,
    fprintf(fh,'%s\n','iter   min    mean     max     STD    reg     rd     tresd');
    temp = [0 min(OPT.DEP.tims) diff(range(OPT.DEP.tims)) max(OPT.DEP.tims) std(OPT.DEP.tims) TST.reg TST.rd TST.tres];
    sprintf(   '%3d %7.1f %7.1f %7.1f  %6.1f %6.1f %7.4f %7.4f',temp(1:8));
    fprintf(fh,'%3d %7.1f %7.1f %7.1f  %6.1f %6.1f %7.4f %7.4f\n',temp(1:8));
    fprintf(fh,'%s\n','iter   min   range  max    std  minAPD rangeAPD maxAPD stdAPD  reg    rd     tres');
end

% disp('iter   min   range  max    STD  minAPD rangeAPD maxAPD stdAPD  reg    rd     tresd');
disp('iter   min   range  max    STD  minAPD maxAPD stdAPD   reg      rd     tresd');

if isempty(diff(range((OPT.DEP.tims)))), dro = NaN; else dro = diff(range((OPT.DEP.tims))); end

temp = [0 min(OPT.DEP.tims) dro max(OPT.DEP.tims) std(OPT.DEP.tims)...
    min(OPT.REP.tims-OPT.DEP.tims) ...
    max(OPT.REP.tims-OPT.DEP.tims) std(OPT.REP.tims-OPT.DEP.tims)...
    TST.reg TST.rd TST.tres];
txt = sprintf(   '%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f Start time: %s',temp,datestr(now,'HH,MM.SS'));
disp(txt);

if INV.doSave,
    fprintf(fh,'%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f',temp);
    fprintf(fh,'time run: %s\n',datestr(now,'HH,MM.SS')) ;
end

RESNOW(1,:) = [temp norm(TST.RES)];
ik          = 2;

%% iterative approach to solving non-linear parameter estimation alternate between dep and rep
iter        = 0;
score       = 1; % start outer loop
startTime   = clock;
depscore    = 1;
repscore    = 1;

if INV.mode == 1, loopstart = 1; nloop = 1; repscore = 0; end   % dep only
if INV.mode == 2, loopstart = 2; nloop = 2; depscore = 0; end   % rep only
if INV.mode == 3, loopstart = 3; nloop = 3; end                 % notch only
if INV.mode == 4, loopstart = 1; nloop = 2; end                 % dep & rep
if INV.mode == 5, loopstart = 1; nloop = 2; end                 % dep & notch & rep ??
if INV.mode == 6, loopstart = 1; nloop = 1; repscore = 0; end   % dep & notch
if INV.mode == 7, loopstart = 7; nloop = 7; end                 % AP amplitude only

totscore    = 0;
TST.rd      = 1;

while iter < 2 || iter < INV.maxiter && (depscore || repscore) &&  TST.rd > INV.MINRD,
    iter = iter + 1;
    for loop = loopstart:nloop,
        switch loop,
            % DEP:
            case 1
                [depscore,OPT.DEP,TST] = optimizeDepRep(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);           %
                if INV.mode == 5 || INV.mode == 6,
                    if INV.doAmplitude,
                        [tmpscore,OPT.AMP,TST] = optimizeAmpl(INV,OPT.AMP,OPT.DEP,OPT.REP,OPT.NOT);     %
                    else
                        [tmpscore,OPT.NOT,TST] = optimizeNotch(INV,OPT.NOT,OPT.DEP,OPT.REP,OPT.AMP);    %
                    end
                end
                opt         = OPT.DEP;
                score       = depscore;
                totscore    = totscore+depscore;
                % REP:
            case 2
                if INV.mode == 8,
                    [tmpscore,OPT.NOT,TST] = optimizeNotch(INV,OPT.NOT,OPT.DEP,OPT.REP,OPT.AMP);
                end
                [repscore,OPT.REP,TST]  = optimizeDepRep(INV,OPT.REP,OPT.DEP,OPT.NOT,OPT.AMP);
                opt                     = OPT.REP;
                score                   = repscore;
                totscore                = totscore+repscore;
                %
            case 3
                [tmpscore,OPT.NOT,TST] = optimizeNotch(INV,OPT.NOT,OPT.DEP,OPT.REP,OPT.AMP);
                opt = OPT.DEP;
                %
            case 7
                [depscore,OPT.DEP,TST]  = optimizeDepRep(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);
                [tmpscore,OPT.AMP,TST]  = optimizeAmpl(INV,OPT.AMP,OPT.DEP,OPT.REP,OPT.NOT);
                [repscore,OPT.REP,TST]  = optimizeDepRep(INV,OPT.REP,OPT.DEP,OPT.NOT,OPT.AMP);
                opt                     = OPT.DEP;
                score                   = depscore;
                totscore                = totscore+depscore;
        end
        
        if isempty(diff(range(opt.tims))), dro = NaN; else dro = diff(range(opt.tims)); end
        
        % no diff-range rep-dep
        temp = [iter min(opt.tims) dro max(opt.tims) std(opt.tims)...
            min(OPT.REP.tims-OPT.DEP.tims) ...
            max(OPT.REP.tims-OPT.DEP.tims) std(OPT.REP.tims-OPT.DEP.tims)...
            TST.reg TST.rd TST.tres];
        
        if score == 0, % ??
            k           = rem(ik,2)+1;
            temp(9:10)  = RESNOW(ik-k,9:10);
            temp(11)    = RESNOW(ik-1,11);
        end
        
        % no diff-range rep-dep
        txt = sprintf('%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.2f %6.3f %6.3f run %s time %s',...
            temp,opt.type,datestr(datenum(clock)-datenum(startTime),'HH,MM.SS'));
        disp(txt);
        
        if INV.doSave,
            fprintf(fh,'%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.2f %6.3f %6.3f',temp);
            fprintf(fh,'run %s time %s\n',opt.type,datestr(datenum(clock)-datenum(startTime),'HH,MM.SS'));
        end
        
        RESNOW(ik,:)    = [temp norm(TST.RES)];
        ik              = ik+1;
    end
end

fprintf('lambda: %f\n',OPT.DEP.lambopt);

%%%% end outer loop; (loop1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TST             = gettres_v(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);
meas.depfinal   = OPT.DEP.tims;
meas.repfinal   = OPT.REP.tims;

if INV.useAmpl,     meas.amplfinal      = OPT.AMP.pot; end
if INV.useNotch,    meas.notchPotfinal  = OPT.NOT.pot; end

meas.rdfinal    = TST.rd;
meas.regfinal   = TST.reg;
meas.tresfinal  = TST.tres;
COR             = corrcoef(TST.PHIA,INV.PHIREF);
meas.corfinal   = COR(2,1);
meas.iterfinal  = iter;
meas.log        = RESNOW;

if INV.doSave, saveasci(outfil,[meas.depfinal meas.repfinal]); fclose(fh); end

%% =======================================================================
function [score,amplopt,TST] = optimizeAmpl(INV,amplopt,depopt,repopt,notchopt)

% prepare compution Marquardt step aimed at improving tims;
% compute Sprime and Sprime*Sprime' based on previous iteration
TST         = gettres_pot(INV,amplopt,depopt,repopt,notchopt);
starttres   = TST.tres;
score       = 0;
lamb        = amplopt.lambopt;

Splus       = getS(INV.T,depopt.tims,repopt.tims,INV.SPECS,notchopt.pot,amplopt.pot+.001,4,INV);
Smin        = getS(INV.T,depopt.tims,repopt.tims,INV.SPECS,notchopt.pot,amplopt.pot-.001,4,INV);
Sprime      = (Splus(:,amplopt.usetimes)-Smin(:,amplopt.usetimes))/2.;
SST         = Sprime*Sprime';

% compute GTGM, with mu setting the weight of the regularization (smoothing based on the surface Laplacian)
GTGM        = INV.ATA.*SST+amplopt.mu^2*INV.ROTRO;

% compute gtres
gtres       = sum(INV.AMA.*(TST.RES*Sprime'))';
righths     = gtres-amplopt.mu^2*INV.ROTRO*amplopt.pot;

% compute next estimate; test residual;
% increase lambda it until starttres < TST.tres
while TST.tres >= starttres, % inner loop
    lamb = lamb*2;
    
    if lamb > INV.lambdamax,  
        break; 
    end
    
    GTGML   = GTGM+lamb^2*eye(size(INV.AMA,2));
    deltau 	= inv(GTGML)*righths;
    
    % 	weight function to limit plateau phase apmlitue to at least 85%
    newpot              = amplopt.pot+deltau;
    newpot(newpot<0)    = 0;
    newpot(newpot>1)    = 1;
    testopt             = amplopt;test.pot=newpot;
    TST                 = gettres_pot(INV,testopt,depopt,repopt,notchopt);
    
    if TST.tres < starttres,
        if (starttres-TST.tres)/starttres >= INV.stopcrit,
            ampl.pot    = newpot;
            score       = 1;
            disp('optimize amplitude')
        end
        
        break;
    end
end  % end of inner loop

amplopt.lambopt = lamb/4;
TST             = gettres_pot(INV,amplopt,depopt,repopt,notchopt);
% TST             = gettres_v(INV,depopt,repopt,notchopt,amplopt);

%% =======================================================================
function [score,notchopt,TST] = optimizeNotch(INV,notchopt,depopt,repopt,amplopt)

% prepare compution Marquardt step aimed at improving tims;
% compute Sprime and Sprime*Sprime' based on previous iteration
TST         = gettres_pot(INV,notchopt,depopt,repopt,amplopt);
starttres   = TST.tres;
score       = 0;
lamb        = notchopt.lambopt;

Splus       = getS(INV.T,depopt.tims,repopt.tims,INV.SPECS,notchopt.pot+0.01,amplopt.pot,4,INV);
Smin        = getS(INV.T,depopt.tims,repopt.tims,INV.SPECS,notchopt.pot-0.01,amplopt.pot,4,INV);
Sprime      = (Splus(:,notchopt.usetimes)-Smin(:,notchopt.usetimes)) / 2.;
SST         = Sprime*Sprime';

% compute GTGM, with mu setting the weight of the regularization
% (smoothing based on the surface Laplacian)
M           = INV.ATA.*SST;
GTGM        = M'*M + + notchopt.mu^2*INV.ROTRO;     % why the + + and not just +
% GTGM       = M + muREG; old
gtres = sum(INV.AMA.*(TST.RES*Sprime'))';

righths     = M'*gtres-muREG*startopt.tims;
righths     = gtres-notchopt.mu^2*INV.ROTRO*notchopt.pot;
testopt     = startopt;

% compute next estimate; test residual;
% increase lambda it until starttres < TST.tres

% compute gtres
righths     = gtres-notchopt.mu^2*INV.ROTRO*notchopt.pot;

% compute next estimate; test residual;
% increase lambda it until starttres < TST.tres

while TST.tres >= starttres, % inner loop
    lamb = lamb*2;
    if lamb > INV.lambdamax,  
        break;
    end
    
    GTGML       = GTGM+lamb^2*eye(size(INV.AMA,2));
    delnotch    = inv(GTGML)*righths;
    
    % 	weight function to limit plateau phase apmlitue to at least 85%
    newpot      = notchopt.pot+delnotch;
    
    % 10 is appoximates the (liniear?) relation between notchpotential
    % 	and the parameter regulating the notch amplitude
    
    newpot(newpot<0)    = 0;
    newpot(newpot>1)    = 1;
    testopt             = notchopt;testopt.pot=newpot;
    TST                 = gettres_pot(INV,testopt,depopt,repopt,amplopt);
    if TST.tres<starttres,
        if (starttres-TST.tres)/starttres >= INV.stopcrit,
            notchopt.pot    = newpot;
            score           = 1;
            disp('optimize notch')
        end
        
        break; 
    end
end  % end of inner loop

notchopt.lambopt    = lamb/4;
TST                 = gettres_pot(INV,notchopt,depopt,repopt,amplopt);
% TST                 = gettres_v(INV,depopt,repopt,notchopt,amplopt);

%% =======================================================================
function [score,opt,TST] = optimizeDepRep(INV,opt,keepopt,notchopt,amplopt)

% prepare compution Marquardt step aimed at improving tims;
% compute Sprime and Sprime*Sprime' based on previous iteration

TST         = gettres_v(INV,opt,keepopt,notchopt,amplopt);
starttres   = TST.tres;
score       = 0;
lamb        = opt.lambopt;

if strcmp(opt.type,'dep'),
    Splus           = getS(INV.T,opt.tims+1,keepopt.tims,INV.SPECS,notchopt.pot,amplopt.pot,INV.mode,INV);
    Smin            = getS(INV.T,opt.tims-1,keepopt.tims,INV.SPECS,notchopt.pot,amplopt.pot,INV.mode,INV);
    startopt        = opt;
elseif strcmp(opt.type,'rep'),
    Splus           = getS(INV.T,keepopt.tims,opt.tims+1,INV.SPECS,notchopt.pot,amplopt.pot,INV.mode,INV);
    Smin            = getS(INV.T,keepopt.tims,opt.tims-1,INV.SPECS,notchopt.pot,amplopt.pot,INV.mode,INV);
    startopt        = opt;
else %if strcmp(opt.type,'apd'),
    Splus           = getS(INV.T,keepopt.tims,opt.tims+1,INV.SPECS,notchopt.pot,amplopt.pot,INV.mode,INV);
    Smin            = getS(INV.T,keepopt.tims,opt.tims-1,INV.SPECS,notchopt.pot,amplopt.pot,INV.mode,INV);
    startopt        = opt;
    startopt.tims   = startopt.tims-keepopt.tims;
end

Sprime  = (Splus-Smin)/2;
SST     = Sprime*Sprime';

% compute GTGM, with mu setting the weight of the regularization
% (smoothing based on the surface Laplacian)
if strcmp(startopt.type,'dep'),
    muREG   = bsxfun(@times,opt.mu^2,INV.ROTRO); % array multiply
    M       = INV.ATA.*SST;
else
    muREG   = bsxfun(@times,opt.mu^2,INV.ROTROREP);
    M       = INV.ATA.*SST;
end

GTGM        = M'*M+muREG;
gtres       = sum(INV.AMA'.*(Sprime * TST.RES'),2); % original PvD
righths     = M'*gtres-muREG*startopt.tims;         % PvD original
testopt     = startopt;

% compute next estimate; test residual;
% increase lambda it until starttres < TST.tres
stuck       = 0; % number of consecutive iterations with small improvements
maxstuck    = 10;

while 1, % TST.tres >= starttres, % inner loop ???
    lamb    = lamb*2;
    
    if lamb > INV.lambdamax, display('lambdamax break'); break; end
    
    GTGM_L          = GTGM+lamb^2*eye(size(INV.AMA,2));     % Orginal AvO: see AvO linsyst.pdf   
    deltau          = GTGM_L\righths;                       % INVERSE; total computing time 0.87 times preceding  
    newtime         = startopt.tims+deltau;    
    testopt.tims    = newtime;
    
    if strcmp(startopt.type,'apd'),
        testopt.tims    = testopt.tims+keepopt.tims;
        TST             = gettres_v(INV,testopt,keepopt,notchopt,amplopt);
    else
        TST             = gettres_v(INV,testopt,keepopt,notchopt,amplopt); 
    end
    
    if TST.tres < starttres, % residual is smaller than in beginning
        if (starttres-TST.tres)/starttres >= INV.stopcrit,
            opt.tims    = testopt.tims;
            score       = 1;
            break
        else
            stuck       = stuck+1;
            if stuck >= maxstuck,
                fprintf('break stuck = %d\n',stuck);
                break;
            end
        end
    else
        stuck = 0; % still optimizing lambda
    end
end  % end of inner loop

opt.lambopt = lamb/4;           % see: AvO linsyst.pdf

%% =======================================================================
function TST = gettres_pot(INV,opt,depopt,repopt,otheropt)

if strcmp(opt.type,'notch'),
    S = getS(INV.T,depopt.tims,repopt.tims,INV.SPECS,opt.pot,otheropt.pot,INV.mode,INV);
else %type == ' ampl',
    S = getS(INV.T,depopt.tims,repopt.tims,INV.SPECS,otheropt.pot,opt.pot,INV.mode,INV);
end

TST.PHIA    = lowpassma(INV.AMA*S(:,opt.usetimes),INV.lpass);
TST.RES     = INV.PHIREF(:,opt.usetimes)-TST.PHIA;
TST.rd      = norm(TST.RES,'fro')/norm(INV.PHIREF(:,opt.usetimes),'fro'); % NOTE: unfiltered rd
TST.reg     = norm(INV.REGOP*opt.pot);%/m2mm;
TST.tres    = sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
TST.S       = S;

%% =======================================================================
function TST = gettres_v(INV,opt,keepopt,notchopt,amplopt)

if strcmp(opt.type,'dep'),
    TST.S = getS(INV.T,opt.tims,keepopt.tims,INV.SPECS,notchopt.pot,amplopt.pot,INV.mode,INV);
else
    TST.S = getS(INV.T,keepopt.tims,opt.tims,INV.SPECS,notchopt.pot,amplopt.pot,INV.mode,INV);
end

TST.PHIA    = lowpassma(INV.AMA*TST.S,INV.lpass);                               % 
TST.RES     = INV.PHIREF-TST.PHIA(1:size(INV.PHIREF,1),1:size(INV.PHIREF,2));   % residue --> difference INV.PHIREF and TST.PHIA
TST.rd      = norm(TST.RES,'fro')/INV.normphi;                                  % NOTE: unfiltered rd
TST.wrd     = sum(rms(TST.RES)./(0.0010+rms(INV.PHIREF)));                      % weighted rd

if strcmp(opt.type,'dep'),
    TST.reg = norm(INV.REGOP*opt.tims); % /m2mm;
else
    TST.reg = norm(INV.REGOPREP*opt.tims); % /m2mm;
end

if INV.useWeighedRd,
    if strcmp(opt.type,'dep'),
        TST.tres = sqrt(TST.wrd^2+(TST.reg*opt.mu)^2);
    else
        TST.tres = sqrt(TST.wrd^2+(TST.reg*opt.mu)^2);
    end
else
    if strcmp(opt.type,'dep'),
        TST.tres = sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
    else
        TST.tres = sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
    end
end

%% =======================================================================
function [REGOP,REGOPREP] = calcREGOP(GEOM,useSurfLapl)

m2mm = 1000;

if useSurfLapl == 1,
    % calculate Laplacian for mesh [GEOM.VER]:
    if max(max(abs(GEOM.VER))) > 1                      % check if vertices of ventricles mesh are in meters or milimeters
        REGOP = surflapl(GEOM.VER/m2mm,GEOM.ITRI,1);
    else
        REGOP = surflapl(GEOM.VER,GEOM.ITRI,1);
    end
    
    REGOPREP = REGOP;
    
else
    disp('inverse 1.distance^2 regularization');
    [a,ORDER]   = graphdist(GEOM.ITRI);
    REGOP       = GEOM.DIST2W;
    
    if max(max(abs(GEOM.VER))) > 1, REGOP = 3*REGOP/m2mm; end % correction if mesh is in meters instead of milimeters
    
    REGOP       = REGOP.^(-1);
    REGOPREP    = REGOP;
    
    if max(max(abs(GEOM.VER))) > 1,
        B = surflapl(GEOM.VER/m2mm,GEOM.ITRI,1);
    else
        B = surflapl(GEOM.VER,GEOM.ITRI,1);
    end
    
    REGOP(GEOM.DIST>10) = 0;
    REGOP(ORDER<=4)     = 0;
    B2                  = diag(-0.3*sum(REGOP));
    REGOP(B~=0)         = B(B~=0);
    REGOP               = REGOP+B2; 
    REGOPREP            = B;
end

%% =======================================================================
function S = getS(T,dep,rep,SPECS,notchpot,scaleAmpl,mode,INV)

if mode == 4,
    SPECS.scaleAmpl = scaleAmpl;
    if size(T,2) < max(rep)+50,
        Tt  = ones(length(dep),1)*(0:max(rep)+50);
        S   = getSmode(Tt,dep,rep,SPECS,4);
        S   = S(:,1:size(T,2));
    else
        S   = getSmode(T,dep,rep,SPECS,4,INV);
    end
else
    SPECS.scaleAmpl = notchpot;
    S               = getSmode(T,dep,rep,SPECS,mode,INV);
end